"""Textual-based terminal UI for editing the pycmplot hits overlay.

This is the "no browser needed" companion to :mod:`pycmplot.editor`
(the Streamlit backend).  It renders in any ANSI terminal, so it
works over vanilla SSH on clusters that don't allow port forwarding
or run headless compute nodes without a graphical display.

Feature parity with the Streamlit backend, adjusted for terminal
ergonomics:

* **DataTable grid** with keyboard-driven navigation (arrow keys /
  ``home`` / ``end`` / ``pgup`` / ``pgdn``).
* **Cell edit modal** — ``enter`` opens an input dialog scoped to the
  focused column.  For ``category`` the dialog remembers existing
  values as an autocomplete-style hint list.  For ``highlight_color``
  a live swatch renders the current value using the terminal's
  truecolor palette (via Rich's ``Style``); invalid values flash a
  red border so typos are caught before save.
* **Save / reload / add / delete row / preview** as footer bindings,
  matching what a spreadsheet user expects.
* **Preview** renders a matplotlib PNG to
  ``<cache_dir>/preview.png`` since most terminals can't render
  images inline; users can ``scp`` or ``sshfs``-open that.  (iTerm2 /
  Kitty could render inline via their graphics protocols but detecting
  them reliably is fiddly and out of scope here.)
* **Quit** prompts for confirmation when there are unsaved changes.

Install with::

    pip install "pycmplot[editor-tui]"

Launch with::

    pycmplot edit --cache_dir ./.pycmplot_cache --tui
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# CLI-side entry point
# ---------------------------------------------------------------------------

def launch(cache_dir: str,
           group: Optional[str] = None) -> None:
    """Boot the Textual editor.  Called from ``pycmplot edit --tui``.

    Fails fast with an actionable message when Textual isn't
    installed, and reuses :func:`pycmplot.editor._resolve_group` so
    the group-picker semantics match the browser backend exactly.
    """
    try:
        import textual  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "The '--tui' backend requires Textual. Install with:\n"
            "    pip install \"pycmplot[editor-tui]\"\n"
            f"(underlying error: {exc})"
        )

    if not Path(cache_dir).exists():
        raise SystemExit(f"cache_dir does not exist: {cache_dir}")

    # Reuse the browser-backend's group resolution — one source of
    # truth for "how do we pick which overlay to edit".
    from pycmplot.editor import _resolve_group
    resolved_group = _resolve_group(cache_dir, group)

    app = HitsEditorApp(cache_dir=str(Path(cache_dir).resolve()),
                        group_key=resolved_group)
    app.run()


# ---------------------------------------------------------------------------
# The App and its screens
# ---------------------------------------------------------------------------

def _build_app_classes():
    """Return the Textual App/Screen classes.

    Defined inside a function so the module imports cleanly even
    without Textual installed — the top-level ``launch()`` can raise
    the helpful install-hint before we ever get here.
    """
    from textual import events, on
    from textual.app import App, ComposeResult
    from textual.binding import Binding
    from textual.containers import Horizontal, Vertical
    from textual.coordinate import Coordinate
    from textual.screen import ModalScreen
    from textual.widgets import (
        Button, DataTable, Footer, Header, Input, Label, Static,
    )

    class _EditModal(ModalScreen[Optional[str]]):
        """Dialog for editing a single cell value.

        Returns the new value (str) on save, or ``None`` on cancel.
        The color-column path renders a live swatch in the modal
        header so you can *see* the value you're picking.  Category
        column shows existing values as a hint line, so autocomplete
        is by copy-paste rather than a full combobox (Textual's
        combobox story is thin as of 8.x and not worth the code).
        """

        BINDINGS = [
            Binding("escape", "cancel", "Cancel", show=True),
            Binding("enter", "submit", "Save", show=True, priority=True),
        ]

        def __init__(self, column: str, current: str,
                     suggestions: Optional[List[str]] = None):
            super().__init__()
            self.column = column
            self.current = current
            self.suggestions = suggestions or []

        def compose(self) -> ComposeResult:
            with Vertical(id="edit-modal"):
                yield Label(f"Edit  [b]{self.column}[/b]", id="edit-title")
                # Live swatch for the colour column — updated on
                # every keystroke via the Input's Changed event.
                if self.column == "highlight_color":
                    yield Static("", id="edit-swatch")
                if self.suggestions:
                    hint = "  ".join(f"[dim]{s}[/dim]"
                                     for s in self.suggestions[:12])
                    yield Static(f"existing: {hint}", id="edit-hint")
                yield Input(value=self.current, id="edit-input")
                with Horizontal(id="edit-actions"):
                    yield Button("Save (⏎)", variant="primary", id="save-btn")
                    yield Button("Cancel (⎋)", id="cancel-btn")

        def on_mount(self) -> None:
            inp = self.query_one("#edit-input", Input)
            inp.focus()
            # Pre-select the current value so typing replaces it —
            # matches the flow spreadsheet users expect (F2 → type
            # → Enter).  ``action_select_all`` exists across Textual
            # 0.60+; if a future version renames it we fall back to
            # a manual cursor-to-end which at least won't corrupt the
            # data.
            try:
                inp.action_select_all()
            except Exception:  # pragma: no cover
                inp.cursor_position = len(inp.value)
            if self.column == "highlight_color":
                self._refresh_swatch(self.current)

        def _refresh_swatch(self, value: str) -> None:
            from matplotlib.colors import is_color_like, to_hex
            v = (value or "").strip()
            if not v or v.lower() in {"auto", "nan", "none", "na"}:
                text = "[dim]auto → falls back to plotter default[/dim]"
            elif is_color_like(v):
                hex_col = to_hex(v)
                # Rich's ``on <hex>`` background makes the terminal
                # emit truecolor SGR — modern terminals render this
                # as an actual coloured chip, older ones fall back
                # to the nearest 256-colour approximation.
                text = f"[white on {hex_col}]   {hex_col}   [/]  valid"
            else:
                text = f"[white on red]   {v}   [/]  [red]invalid[/red]"
            self.query_one("#edit-swatch", Static).update(text)

        @on(Input.Changed, "#edit-input")
        def _on_change(self, event: Input.Changed) -> None:
            if self.column == "highlight_color":
                self._refresh_swatch(event.value)

        def action_submit(self) -> None:
            self.dismiss(self.query_one("#edit-input", Input).value)

        def action_cancel(self) -> None:
            self.dismiss(None)

        @on(Button.Pressed, "#save-btn")
        def _save(self) -> None:
            self.action_submit()

        @on(Button.Pressed, "#cancel-btn")
        def _cancel(self) -> None:
            self.action_cancel()


    class _ConfirmModal(ModalScreen[bool]):
        """Yes/No modal used for quit-with-unsaved and delete-row."""

        BINDINGS = [
            Binding("escape", "cancel", "No", show=True),
            Binding("y", "confirm", "Yes", show=True),
            Binding("n", "cancel", "No", show=False),
        ]

        def __init__(self, message: str):
            super().__init__()
            self.message = message

        def compose(self) -> ComposeResult:
            with Vertical(id="confirm-modal"):
                yield Label(self.message, id="confirm-msg")
                with Horizontal(id="confirm-actions"):
                    yield Button("Yes (y)", variant="warning", id="yes")
                    yield Button("No (n)", id="no")

        def action_confirm(self) -> None:
            self.dismiss(True)

        def action_cancel(self) -> None:
            self.dismiss(False)

        @on(Button.Pressed, "#yes")
        def _yes(self) -> None: self.action_confirm()

        @on(Button.Pressed, "#no")
        def _no(self) -> None: self.action_cancel()


    class _FilterModal(ModalScreen[Optional[str]]):
        """Prompt for a pandas ``df.query()`` expression.

        The empty string is a valid submission — it clears the filter,
        matching the CLI's semantics for "no --where".
        """

        BINDINGS = [
            Binding("escape", "cancel", "Cancel", show=True),
            Binding("enter", "submit", "Apply", show=True, priority=True),
        ]

        def __init__(self, current: str = ""):
            super().__init__()
            self.current = current

        def compose(self) -> ComposeResult:
            with Vertical(id="edit-modal"):
                yield Label("Filter  [dim](pandas df.query grammar)[/dim]",
                            id="edit-title")
                yield Static(
                    "[dim]Examples:  P < 5e-8   |   "
                    "category == \"significant\"   |   "
                    "CHR == \"6\" and BP.between(28e6, 34e6)[/dim]",
                    id="edit-hint",
                )
                yield Input(value=self.current,
                            placeholder="empty submission clears filter",
                            id="edit-input")
                with Horizontal(id="edit-actions"):
                    yield Button("Apply (⏎)", variant="primary", id="save-btn")
                    yield Button("Cancel (⎋)", id="cancel-btn")

        def on_mount(self) -> None:
            inp = self.query_one("#edit-input", Input)
            inp.focus()
            try:
                inp.action_select_all()
            except Exception:  # pragma: no cover
                inp.cursor_position = len(inp.value)

        def action_submit(self) -> None:
            self.dismiss(self.query_one("#edit-input", Input).value)

        def action_cancel(self) -> None:
            self.dismiss(None)

        @on(Button.Pressed, "#save-btn")
        def _save(self) -> None: self.action_submit()

        @on(Button.Pressed, "#cancel-btn")
        def _cancel(self) -> None: self.action_cancel()


    class _BatchEditModal(ModalScreen[Optional[Tuple[str, str]]]):
        """Pick a column (color/category) and a value; apply to N rows.

        Returns ``(column_name, value)`` on save, or ``None`` on
        cancel.  The color column path renders the same live swatch
        as :class:`_EditModal` so users see the resolved chip before
        committing.
        """

        BINDINGS = [
            Binding("escape", "cancel", "Cancel", show=True),
            Binding("enter", "submit", "Apply", show=True, priority=True),
            # ``tab`` switches column focus; we hint it so first-time
            # users know they can flip between highlight_color and
            # category without reaching for the mouse.
            Binding("tab", "focus_next", "Next field", show=False),
        ]

        def __init__(self, n_rows: int,
                     category_hints: Optional[List[str]] = None):
            super().__init__()
            self.n_rows = n_rows
            self.category_hints = category_hints or []
            self._column = "highlight_color"

        def compose(self) -> ComposeResult:
            with Vertical(id="edit-modal"):
                yield Label(
                    f"Batch edit  [b]{self.n_rows}[/] row(s)",
                    id="edit-title",
                )
                with Horizontal(id="col-toggle"):
                    yield Button(
                        "highlight_color", id="pick-color",
                        variant="primary",
                    )
                    yield Button("category", id="pick-category")
                yield Static("", id="edit-swatch")
                if self.category_hints:
                    hint = "  ".join(f"[dim]{s}[/dim]"
                                     for s in self.category_hints[:12])
                    yield Static(f"existing: {hint}", id="edit-hint")
                yield Input(value="", id="edit-input",
                            placeholder="value to write into every selected row")
                with Horizontal(id="edit-actions"):
                    yield Button(
                        "Apply (⏎)", variant="primary", id="save-btn",
                    )
                    yield Button("Cancel (⎋)", id="cancel-btn")

        def on_mount(self) -> None:
            self.query_one("#edit-input", Input).focus()
            self._refresh_swatch("")

        def _refresh_swatch(self, value: str) -> None:
            widget = self.query_one("#edit-swatch", Static)
            if self._column != "highlight_color":
                widget.update("")
                return
            from matplotlib.colors import is_color_like, to_hex
            v = (value or "").strip()
            if not v or v.lower() in {"auto", "nan", "none", "na"}:
                text = "[dim]auto → falls back to plotter default[/dim]"
            elif is_color_like(v):
                hex_col = to_hex(v)
                text = f"[white on {hex_col}]   {hex_col}   [/]  valid"
            else:
                text = f"[white on red]   {v}   [/]  [red]invalid[/red]"
            widget.update(text)

        @on(Input.Changed, "#edit-input")
        def _on_change(self, event: Input.Changed) -> None:
            self._refresh_swatch(event.value)

        @on(Button.Pressed, "#pick-color")
        def _pick_color(self) -> None:
            self._column = "highlight_color"
            self.query_one("#pick-color", Button).variant = "primary"
            self.query_one("#pick-category", Button).variant = "default"
            self._refresh_swatch(self.query_one("#edit-input", Input).value)

        @on(Button.Pressed, "#pick-category")
        def _pick_category(self) -> None:
            self._column = "category"
            self.query_one("#pick-color", Button).variant = "default"
            self.query_one("#pick-category", Button).variant = "primary"
            self._refresh_swatch("")  # category has no swatch

        def action_submit(self) -> None:
            val = self.query_one("#edit-input", Input).value.strip()
            if not val:
                # Empty value is a no-op — dismiss without committing
                # rather than clobbering every selected row with "".
                self.dismiss(None)
                return
            if self._column == "highlight_color":
                from matplotlib.colors import is_color_like
                if val.lower() != "auto" and not is_color_like(val):
                    # Refuse the write; let the user see the red
                    # swatch and try again.
                    return
            self.dismiss((self._column, val))

        def action_cancel(self) -> None:
            self.dismiss(None)

        @on(Button.Pressed, "#save-btn")
        def _save(self) -> None: self.action_submit()

        @on(Button.Pressed, "#cancel-btn")
        def _cancel(self) -> None: self.action_cancel()


    class HitsEditorApp(App[None]):
        """Main TUI application.

        Holds the working dataframe in memory (``self.df``) and only
        writes to disk when the user hits Save.  The dataframe is the
        single source of truth; the DataTable widget is refreshed
        from it after every mutation.
        """

        # Some minimal CSS keeps the modals from spanning the whole
        # terminal.  Textual's CSS is a subset of web CSS — no external
        # file so the app is single-file self-contained.
        CSS = """
        #edit-modal, #confirm-modal {
            align: center middle;
            width: 60;
            height: auto;
            border: round $accent;
            padding: 1 2;
            background: $panel;
        }
        #edit-title, #confirm-msg { margin-bottom: 1; }
        #edit-swatch, #edit-hint { margin-bottom: 1; }
        #edit-actions, #confirm-actions {
            height: 3;
            margin-top: 1;
            align: center middle;
        }
        /* Batch-edit column toggle: two buttons side-by-side above
           the value input. */
        #col-toggle {
            height: 3;
            margin-bottom: 1;
            align: center middle;
        }
        Button { margin: 0 1; }
        #status { dock: bottom; height: 1; background: $primary-darken-2; }
        """

        BINDINGS = [
            Binding("ctrl+s", "save", "Save", show=True),
            Binding("ctrl+r", "reload", "Reload", show=True),
            Binding("ctrl+n", "add_row", "Add row", show=True),
            Binding("ctrl+d", "delete_row", "Delete row", show=True),
            Binding("ctrl+p", "preview", "Preview PNG", show=True),
            # Filter + multi-select + batch-edit — the "hundreds of
            # loci" flow.  ``/`` opens a filter dialog (pandas
            # ``df.query()`` grammar, same as the CLI's ``--where``).
            # ``space`` toggles selection on the focused row.
            # ``ctrl+a`` selects every row currently visible.
            # ``ctrl+e`` opens a batch-edit modal that applies one
            # colour or category value to every selected row.
            # ``escape`` clears filter + selection so the user gets
            # the full view back with one key.
            Binding("slash", "filter", "Filter", show=True),
            Binding("space", "toggle_select", "Select", show=True),
            Binding("ctrl+a", "select_all", "Select all", show=True),
            Binding("ctrl+e", "batch_edit", "Batch edit", show=True),
            Binding("escape", "clear_filter_and_selection", "Clear",
                    show=False),
            # DataTable already binds ``enter`` to its own
            # ``CellSelected`` action, so a priority app-level binding
            # here would collide.  We use ``f2`` (Excel convention)
            # as the primary shortcut and additionally hook
            # ``DataTable.CellSelected`` below so ``enter`` still opens
            # the edit modal — the user gets both keys without any
            # binding-precedence surprises.
            Binding("f2", "edit_cell", "Edit", show=True),
            Binding("ctrl+q", "quit", "Quit", show=True),
        ]

        def __init__(self, cache_dir: str, group_key: str):
            super().__init__()
            self.cache_dir = cache_dir
            self.group_key = group_key
            self.df = None  # type: ignore[assignment]
            self.meta: Optional[dict] = None
            self.dirty: bool = False
            self._status: Optional[Static] = None
            # Filter + selection state.
            #
            # ``filter_expr`` is a pandas df.query() string; ``None``
            # means show everything.  ``selected`` holds the *original*
            # dataframe indices of rows the user has multi-selected —
            # keeping them tied to df indices rather than grid rows
            # means the selection survives filtering and refresh.
            # ``_display_indices`` maps ``grid_row -> df_index`` so
            # cursor-based actions know which underlying row to
            # mutate.
            self.filter_expr: Optional[str] = None
            self.selected: set = set()
            self._display_indices: List[int] = []

        # ---- compose / layout -------------------------------------
        def compose(self) -> ComposeResult:
            yield Header(show_clock=False)
            yield DataTable(id="grid", zebra_stripes=True, cursor_type="cell")
            yield Static("", id="status")
            yield Footer()

        # ---- lifecycle --------------------------------------------
        def on_mount(self) -> None:
            self.title = "pycmplot hits editor (TUI)"
            self.sub_title = (
                f"{Path(self.cache_dir).name}/annotations/"
                f"hits.{self.group_key}.tsv"
            )
            self._status = self.query_one("#status", Static)
            self._reload_from_disk()

        # ---- data plumbing ---------------------------------------
        def _reload_from_disk(self) -> None:
            from pycmplot.cache import read_hits_overlay
            df, meta = read_hits_overlay(self.cache_dir, self.group_key)
            if df is None:
                self._notify(
                    f"Could not read hits.{self.group_key}.tsv "
                    "(re-run pycmplot --cache to regenerate)."
                )
                return
            self.df = df.reset_index(drop=True)
            self.meta = meta
            self.dirty = False
            # A fresh dataframe invalidates any filter that referenced
            # by-index rows, and blows away selections that pointed at
            # the previous index space.
            self.filter_expr = None
            self.selected.clear()
            self._refresh_grid()
            self._notify(f"Loaded {len(self.df)} row(s).")

        def _get_visible_df(self):
            """Return the currently-visible (filtered) view of ``self.df``.

            The returned frame keeps the original df's index, which
            is what ``_display_indices`` is built from — never
            ``reset_index`` here or the selection/mutation code loses
            its mapping.
            """
            if self.df is None or self.df.empty:
                return self.df
            if not self.filter_expr:
                return self.df
            try:
                return self.df.query(self.filter_expr, engine="python")
            except Exception as exc:
                self._notify(f"Filter error: {exc}")
                return self.df

        def _refresh_grid(self) -> None:
            grid = self.query_one("#grid", DataTable)
            grid.clear(columns=True)
            if self.df is None or self.df.empty:
                grid.add_column("(empty)", width=20)
                self._display_indices = []
                return
            visible = self._get_visible_df()
            # First column is the selection marker.  It's not part of
            # the dataframe — we render it from ``self.selected`` so
            # the marker survives sort/filter/refresh.
            grid.add_column("✓", key="__sel__", width=3)
            for col in self.df.columns:
                grid.add_column(str(col), key=str(col))
            self._display_indices = list(visible.index)
            for df_idx, row in visible.iterrows():
                sel = "[b green]●[/]" if df_idx in self.selected else " "
                cells = [sel] + [
                    self._render_cell(col, row[col])
                    for col in self.df.columns
                ]
                grid.add_row(*cells, key=str(df_idx))

        @staticmethod
        def _render_cell(col: str, val: Any) -> str:
            """Render a cell value.

            Colour cells get a coloured chip prepended so the user
            sees the mapping without having to open the edit modal;
            other cells are rendered as plain strings.
            """
            s = "" if val is None else str(val)
            if col == "highlight_color":
                from matplotlib.colors import is_color_like, to_hex
                v = s.strip()
                if not v or v.lower() in {"auto", "nan", "none"}:
                    return f"[dim]{s or 'auto'}[/dim]"
                if is_color_like(v):
                    return f"[on {to_hex(v)}]  [/] {s}"
                return f"[white on red] ! [/] {s}"
            return s

        # ---- events -----------------------------------------------
        @on(DataTable.CellSelected)
        def _on_cell_selected(self, event: DataTable.CellSelected) -> None:
            # ``DataTable`` fires this on ``enter`` (and double-click);
            # forwarding to ``action_edit_cell`` gives users the
            # familiar spreadsheet flow without stealing the binding
            # from DataTable's own key handling.
            self.action_edit_cell()

        # ---- cursor -> (df_index, col_name) translation ----------
        def _cursor_target(self) -> Optional[Tuple[int, str]]:
            """Return the ``(df_index, column_name)`` under the cursor.

            Accounts for the leading selection-marker column (grid col
            0) and for the fact that grid rows are indexed into the
            filtered view via ``_display_indices``.  Returns ``None``
            when the cursor is on the marker column or outside the
            data area.
            """
            grid = self.query_one("#grid", DataTable)
            if grid.cursor_coordinate is None or self.df is None:
                return None
            g_row = grid.cursor_coordinate.row
            g_col = grid.cursor_coordinate.column
            if g_row < 0 or g_row >= len(self._display_indices):
                return None
            if g_col == 0:  # marker column
                return None
            data_col = g_col - 1
            if data_col >= len(self.df.columns):
                return None
            return self._display_indices[g_row], str(self.df.columns[data_col])

        # ---- actions ----------------------------------------------
        def action_edit_cell(self) -> None:
            target = self._cursor_target()
            if target is None:
                return
            df_idx, col_name = target
            current = "" if self.df.at[df_idx, col_name] is None \
                else str(self.df.at[df_idx, col_name])

            suggestions: List[str] = []
            if col_name == "category":
                suggestions = sorted({str(v) for v in self.df[col_name].dropna()
                                      if str(v).strip()})

            def _on_close(new_val: Optional[str]) -> None:
                if new_val is None:
                    return
                self.df.at[df_idx, col_name] = new_val
                self.dirty = True
                self._refresh_grid()
                self._notify(f"Set {col_name}[{df_idx}] = {new_val!r}  (unsaved)")

            self.push_screen(_EditModal(col_name, current, suggestions),
                             _on_close)

        # ---- new: filter -----------------------------------------
        def action_filter(self) -> None:
            """Open a modal Input to enter a pandas ``df.query()`` expression.

            The expression is stored on ``self.filter_expr`` and
            applied by ``_get_visible_df`` on every grid refresh, so
            selection state (indexed on the *original* df) survives.
            Submitting an empty expression clears the filter.
            """
            def _on_close(new_val: Optional[str]) -> None:
                if new_val is None:
                    return
                expr = new_val.strip()
                self.filter_expr = expr or None
                self._refresh_grid()
                visible = self._get_visible_df()
                n_vis = len(visible) if visible is not None else 0
                total = 0 if self.df is None else len(self.df)
                if self.filter_expr:
                    self._notify(
                        f"Filter: {self.filter_expr!r}  ({n_vis}/{total})"
                    )
                else:
                    self._notify(f"Filter cleared  ({total} row(s))")

            self.push_screen(
                _FilterModal(current=self.filter_expr or ""),
                _on_close,
            )

        # ---- new: multi-select -----------------------------------
        def action_toggle_select(self) -> None:
            grid = self.query_one("#grid", DataTable)
            if grid.cursor_coordinate is None or self.df is None:
                return
            g_row = grid.cursor_coordinate.row
            if g_row < 0 or g_row >= len(self._display_indices):
                return
            df_idx = self._display_indices[g_row]
            if df_idx in self.selected:
                self.selected.discard(df_idx)
            else:
                self.selected.add(df_idx)
            self._refresh_grid()
            grid.move_cursor(row=g_row, column=grid.cursor_coordinate.column)
            self._notify(
                f"{len(self.selected)} row(s) selected."
            )

        def action_select_all(self) -> None:
            """Select every row currently visible under the filter."""
            if self.df is None:
                return
            visible = self._get_visible_df()
            if visible is None or visible.empty:
                return
            self.selected.update(visible.index)
            self._refresh_grid()
            self._notify(f"{len(self.selected)} row(s) selected.")

        def action_clear_filter_and_selection(self) -> None:
            had_filter = self.filter_expr is not None
            had_sel = bool(self.selected)
            self.filter_expr = None
            self.selected.clear()
            if had_filter or had_sel:
                self._refresh_grid()
                self._notify("Filter + selection cleared.")

        # ---- new: batch edit -------------------------------------
        def action_batch_edit(self) -> None:
            """Apply one colour or category value to every selected row.

            Falls back to "the current filtered view" when no rows
            are explicitly selected — that way ``/`` + ``ctrl+e`` is a
            two-step batch flow, not three.
            """
            if self.df is None:
                return
            targets = list(self.selected)
            if not targets:
                visible = self._get_visible_df()
                if visible is None or visible.empty:
                    self._notify("Nothing to edit.")
                    return
                targets = list(visible.index)
                _source = "current filtered view"
            else:
                _source = "current selection"

            # Category suggestions for the modal's hint line.
            cats = sorted({str(v) for v in self.df.get("category", []).dropna()
                           if str(v).strip()})

            def _on_close(result: Optional[Tuple[str, str]]) -> None:
                if result is None:
                    return
                col, val = result
                self.df.loc[targets, col] = val
                self.dirty = True
                self._refresh_grid()
                self._notify(
                    f"Batch-set {col} = {val!r} on {len(targets)} row(s) "
                    f"({_source}).  Unsaved."
                )

            self.push_screen(
                _BatchEditModal(n_rows=len(targets), category_hints=cats),
                _on_close,
            )

        def action_save(self) -> None:
            if self.df is None:
                self._notify("Nothing to save.")
                return
            if not self.meta or "auto_key" not in (self.meta or {}):
                self._notify(
                    "Missing auto_key in sidecar; refusing save "
                    "(re-run pycmplot --cache to regenerate meta)."
                )
                return
            from pycmplot.cache import (
                AUTO_TAG, SOURCE_COL, write_hits_overlay,
            )
            try:
                auto_df = self.df[
                    self.df[SOURCE_COL].astype(str) == AUTO_TAG
                ].copy()
                written = write_hits_overlay(
                    self.cache_dir, auto_df,
                    auto_key=self.meta["auto_key"],
                    group_key=self.group_key,
                    preserve_user_from=self.df,
                )
                self.dirty = False
                self._notify(f"Saved: {written}")
            except Exception as exc:  # pragma: no cover — real error path
                self._notify(f"Save failed: {exc}")

        def action_reload(self) -> None:
            if self.dirty:
                def _cb(yes: bool) -> None:
                    if yes:
                        self._reload_from_disk()
                self.push_screen(
                    _ConfirmModal("Discard unsaved changes and reload?"),
                    _cb,
                )
            else:
                self._reload_from_disk()

        def action_add_row(self) -> None:
            if self.df is None:
                return
            import pandas as pd
            from pycmplot.annotation import (
                CATEGORY_COL, CATEGORY_DEFAULT,
                HIGHLIGHT_COLOR_AUTO, HIGHLIGHT_COLOR_COL,
            )
            from pycmplot.cache import SOURCE_COL, USER_TAG

            blank = {c: "" for c in self.df.columns}
            if SOURCE_COL in blank:
                blank[SOURCE_COL] = USER_TAG
            if HIGHLIGHT_COLOR_COL in blank:
                blank[HIGHLIGHT_COLOR_COL] = HIGHLIGHT_COLOR_AUTO
            if CATEGORY_COL in blank:
                blank[CATEGORY_COL] = CATEGORY_DEFAULT
            self.df = pd.concat(
                [self.df, pd.DataFrame([blank])], ignore_index=True,
            )
            self.dirty = True
            self._refresh_grid()
            self.query_one("#grid", DataTable).move_cursor(
                row=len(self.df) - 1, column=0,
            )
            self._notify(f"Added user row (index {len(self.df) - 1}).")

        def action_delete_row(self) -> None:
            grid = self.query_one("#grid", DataTable)
            if grid.cursor_coordinate is None or self.df is None:
                return
            g_row = grid.cursor_coordinate.row
            if g_row < 0 or g_row >= len(self._display_indices):
                return
            df_idx = self._display_indices[g_row]

            def _cb(yes: bool) -> None:
                if not yes:
                    return
                self.df = self.df.drop(df_idx).reset_index(drop=True)
                # Deleting a row breaks the selection index space —
                # simplest is to blow selection away.  Filter is fine
                # because df.query() operates on the new frame.
                self.selected.clear()
                self.dirty = True
                self._refresh_grid()
                self._notify(f"Deleted 1 row. (unsaved)")

            self.push_screen(
                _ConfirmModal(f"Delete this row?  Not undoable."),
                _cb,
            )

        def action_preview(self) -> None:
            """Render a linear Manhattan PNG using in-memory overlay."""
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from pycmplot.cache import TrackCache
            from pycmplot import __version__ as _pcm_version

            if self.df is None:
                self._notify("Nothing to preview.")
                return

            tc = TrackCache(self.cache_dir, _pcm_version)
            tracks = getattr(tc, "_metadata", {}).get("tracks", {}) or {}
            by_label: Dict[str, Tuple[str, dict]] = {}
            for ckey, entry in tracks.items():
                lbl = entry.get("label")
                if not lbl:
                    continue
                prev = by_label.get(lbl)
                if prev is None or entry.get("written", 0) > prev[1].get("written", 0):
                    by_label[lbl] = (ckey, entry)
            if not by_label:
                self._notify("No cached tracks — run pycmplot --cache first.")
                return
            # ``plot_linear`` expects ``sumstats_loaded`` to be a
            # dict of ``label -> [df, n_chroms]`` (list, indexed
            # positionally at line 1699 of ``pycmplot/plotting/linear.py``).
            # A previous version of this preview passed raw DataFrames
            # and hit ``KeyError(0)`` — restoring the tuple shape here.
            dfs = {}
            for lbl, (ck, _e) in by_label.items():
                got = tc.get(lbl, ck)
                if got is not None and getattr(got, "df", None) is not None:
                    n_chroms = int(getattr(got, "extras", {}).get(
                        "n_chroms", got.df["CHR"].nunique() if "CHR" in got.df.columns else 22,
                    ))
                    dfs[lbl] = [got.df, n_chroms]
            if not dfs:
                self._notify("Cache metadata present but parquet files unreadable.")
                return

            fig, _ = plt.subplots(figsize=(12, 3 * max(1, len(dfs))))
            try:
                from pycmplot.plotting.linear import linear
                linear(
                    sumstats_loaded=dfs, logp=True, highlight=True,
                    hits_table=self.df, annotate=False, ax=fig,
                )
                out_path = Path(self.cache_dir) / "preview.png"
                fig.savefig(out_path, dpi=140, bbox_inches="tight")
                plt.close(fig)
                self._notify(f"Preview saved: {out_path}")
            except Exception as exc:  # pragma: no cover
                plt.close(fig)
                self._notify(f"Preview failed: {exc}")

        def action_quit(self) -> None:
            if not self.dirty:
                self.exit()
                return

            def _cb(yes: bool) -> None:
                if yes:
                    self.exit()

            self.push_screen(
                _ConfirmModal("Unsaved changes.  Quit anyway?"),
                _cb,
            )

        # ---- helpers ---------------------------------------------
        def _notify(self, msg: str) -> None:
            if self._status is not None:
                self._status.update(msg)
            logger.info("pycmplot editor_tui: %s", msg)

    return HitsEditorApp


# ``launch()`` needs the App class; lazy-build it so that importing
# this module never touches Textual (which may be absent).
def HitsEditorApp(*args, **kwargs):
    return _build_app_classes()(*args, **kwargs)
