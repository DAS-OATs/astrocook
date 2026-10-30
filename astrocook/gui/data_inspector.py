"""The data inspector: a table of the pixels around a point of the main
plot, opened from its context menu, with a sliding window of rows.

Moved out of ``pyside_plot.py`` in phase 5 (split of the GUI).
"""
import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (QAbstractItemView, QDialog, QDialogButtonBox, QHeaderView,
                               QTableWidget, QTableWidgetItem, QVBoxLayout)
from .theme import SELECTION


class DataInspectorDialog(QDialog):
    """
    Non-blocking popup table with Sliding Window scrolling.
    - Infinite scrolling in both directions.
    - Constant memory usage (prunes rows as you scroll).
    - Persistent highlighting of the original center row.
    """
    def __init__(self, session, center_idx, window=20, parent=None):
        super().__init__(parent)
        
        self.session = session
        self.spec = session.spec
        self.center_idx = center_idx  # [MEMORY] Store original target
        
        # Sliding Window Config
        self.chunk_size = 20     # Rows to add per load
        self.max_rows = 200      # Max rows to keep in memory before pruning
        
        # Track current range
        self.current_start = max(0, center_idx - window)
        self.current_end = min(len(self.spec.x), center_idx + window + 1)
        
        # View Setup
        x_val = self.spec.x.value[center_idx]
        x_unit = str(self.spec.x.unit)
        self.setWindowTitle(f"Data Inspector (x = {x_val:.5f} {x_unit})")
        self.resize(700, 450)
        
        layout = QVBoxLayout(self)
        border_width = 10 
        layout.setContentsMargins(border_width, border_width, border_width, border_width)

        # 1. Setup Table
        self.table = QTableWidget()
        
        # Define Columns
        self.cols = ['x', 'xmin', 'xmax', 'y', 'dy', 'cont', 'model']

        # Define Display Labels
        self.display_labels = {
            'x': 'λ',
            'xmin': 'λmin',
            'xmax': 'λmax',
            'y': 'F',
            'dy': 'dF',
            'cont': 'cont',
            'model': 'model'
        }
        
        if hasattr(self.spec, '_data') and self.spec._data.aux_cols:
            for c in self.spec._data.aux_cols:
                if c not in self.cols: 
                    self.cols.append(c)
                    self.display_labels[c] = c # Default label for aux columns
                
        self.table.setColumnCount(len(self.cols) + 1)
        
        # Use the display_labels for the horizontal header
        header_labels = ['Index'] + [self.display_labels.get(c, c) for c in self.cols]
        self.table.setHorizontalHeaderLabels(header_labels)
        
        # Styling
        self.table.verticalHeader().setVisible(False)
        header_font = self.table.horizontalHeader().font()
        header_font.setBold(True)
        self.table.horizontalHeader().setFont(header_font)
        self.table.setAlternatingRowColors(True)
        
        # 2. Initial Population
        self._populate_range(self.current_start, self.current_end, prepend=False)
        
        # 3. Scroll & Layout
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        layout.addWidget(self.table)
        
        btn_box = QDialogButtonBox(QDialogButtonBox.Close)
        btn_box.rejected.connect(self.close)
        layout.addWidget(btn_box)

        # 4. Scroll to Center initially
        center_row = center_idx - self.current_start
        if 0 <= center_row < self.table.rowCount():
            item = self.table.item(center_row, 0)
            self.table.scrollToItem(item, QAbstractItemView.PositionAtCenter)
            
        # 5. Connect Scroll Signal
        self.table.verticalScrollBar().valueChanged.connect(self._check_scroll)

    def refresh(self, new_session):
        """
        Updates the table data and columns if the session state changes 
        (e.g., a column is added or masked), while preserving the scroll window.
        """
        if not new_session or not new_session.spec: 
            return
            
        self.session = new_session
        self.spec = new_session.spec
        
        # 1. Save current scroll position to prevent jumping
        vbar = self.table.verticalScrollBar()
        scroll_val = vbar.value()
        
        # 2. Re-evaluate columns in case one was added or removed
        self.cols = ['x', 'xmin', 'xmax', 'y', 'dy', 'cont', 'model']
        self.display_labels = {
            'x': 'λ', 'xmin': 'λmin', 'xmax': 'λmax', 
            'y': 'F', 'dy': 'dF', 'cont': 'cont', 'model': 'model'
        }
        
        if hasattr(self.spec, '_data') and self.spec._data.aux_cols:
            for c in self.spec._data.aux_cols:
                if c not in self.cols: 
                    self.cols.append(c)
                    self.display_labels[c] = c # Default label for aux columns
                    
        # 3. Block signals to prevent scroll events from firing during rebuild
        self.table.blockSignals(True)
        vbar.blockSignals(True)
        
        # 4. Rebuild table headers and clear old rows
        self.table.setColumnCount(len(self.cols) + 1)
        header_labels = ['Index'] + [self.display_labels.get(c, c) for c in self.cols]
        self.table.setHorizontalHeaderLabels(header_labels)
        self.table.setRowCount(0) 
        
        self.table.blockSignals(False)
        
        # 5. Re-populate the current sliding window with the new data
        self._populate_range(self.current_start, self.current_end, prepend=False)
        
        # 6. Restore scroll position
        vbar.setValue(scroll_val)
        vbar.blockSignals(False)

    def _populate_range(self, start_idx, end_idx, prepend=False):
        """
        Inserts rows. Checks every row against self.center_idx to apply highlight.
        """
        indices = range(start_idx, end_idx)
        if not indices: return

        self.table.setSortingEnabled(False)
        self.table.blockSignals(True) # Optimization
        
        for i, data_idx in enumerate(indices):
            if prepend:
                row_pos = i 
            else:
                row_pos = self.table.rowCount()
                
            self.table.insertRow(row_pos)
            
            # --- Check Highlight ---
            is_highlight = (data_idx == self.center_idx)
            
            # Index Column
            item_idx = QTableWidgetItem(str(data_idx))
            item_idx.setTextAlignment(Qt.AlignCenter)
            if is_highlight:
                self._apply_highlight(item_idx)
            self.table.setItem(row_pos, 0, item_idx)

            # Data Columns
            for col_i, col_name in enumerate(self.cols):
                val_str = "N/A"
                try:
                    if hasattr(self.spec, col_name):
                        col_obj = getattr(self.spec, col_name)
                        val = col_obj.value[data_idx] if col_obj is not None else np.nan
                    else:
                        col_obj = self.spec.get_column(col_name)
                        val = col_obj.value[data_idx] if col_obj is not None else np.nan
                        
                    if isinstance(val, (float, np.floating)):
                        val_str = f"{val:.5e}"
                    else:
                        val_str = str(val)
                except:
                    pass

                item = QTableWidgetItem(val_str)
                item.setTextAlignment(Qt.AlignCenter)
                if is_highlight:
                    self._apply_highlight(item)
                self.table.setItem(row_pos, col_i + 1, item)
        
        self.table.blockSignals(False)

    def _apply_highlight(self, item):
        """Helper to apply the persistent orange highlight."""
        item.setBackground(QColor(SELECTION))
        font = item.font()
        font.setBold(True)
        item.setFont(font)

    def _check_scroll(self, value):
        """
        Sliding Window Logic: 
        - Adds rows to one end.
        - Prunes from the other.
        - Counter-adjusts scrollbar to maintain visual stability ('bounce' effect).
        """
        vbar = self.table.verticalScrollBar()
        
        # Define Thresholds (trigger before hitting the absolute edge for smoothness)
        # using a small buffer (e.g. 10 pixels) or just check min/max
        is_top = value <= vbar.minimum()
        is_bottom = value >= vbar.maximum()
        
        if not (is_top or is_bottom):
            return

        # Optimization: Prevent recursive calls while modifying table
        self.table.blockSignals(True)
        vbar.blockSignals(True)
        
        try:
            if is_bottom:
                # 1. Append New Chunk to Bottom
                new_start = self.current_end
                new_end = min(len(self.spec.x), self.current_end + self.chunk_size)
                
                if new_end > new_start:
                    self._populate_range(new_start, new_end, prepend=False)
                    self.current_end = new_end
                    
                    # 2. Prune Top (Sliding Window)
                    if self.table.rowCount() > self.max_rows:
                        rows_to_remove = self.chunk_size
                        
                        # [FIX] Measure height BEFORE removing to be exact
                        # (Summing height of rows 0..N ensures accuracy if row heights vary)
                        pixels_removed = 0
                        for r in range(rows_to_remove):
                            pixels_removed += self.table.rowHeight(r)
                            
                        # Remove rows from top
                        for _ in range(rows_to_remove):
                            self.table.removeRow(0)
                        self.current_start += rows_to_remove
                        
                        # [FIX] Counter-adjust Scrollbar (The "Bounce Back Up")
                        # The data we were looking at moved UP, so we must scroll UP.
                        vbar.setValue(vbar.value() - pixels_removed)

            elif is_top:
                # 1. Prepend New Chunk to Top
                new_end = self.current_start
                new_start = max(0, self.current_start - self.chunk_size)
                
                if new_end > new_start:
                    self._populate_range(new_start, new_end, prepend=True)
                    self.current_start = new_start
                    
                    # 2. Counter-adjust Scrollbar (The "Bounce Down")
                    # We added rows to the top, pushing our data DOWN.
                    # We must scroll DOWN to keep looking at the same data.
                    rows_added = new_end - new_start
                    pixels_added = 0
                    for r in range(rows_added):
                        pixels_added += self.table.rowHeight(r)
                        
                    vbar.setValue(vbar.value() + pixels_added)
                    
                    # 3. Prune Bottom (Sliding Window)
                    if self.table.rowCount() > self.max_rows:
                        rows_to_remove = self.chunk_size
                        for _ in range(rows_to_remove):
                            self.table.removeRow(self.table.rowCount() - 1)
                        self.current_end -= rows_to_remove
                        
        finally:
            self.table.blockSignals(False)
            vbar.blockSignals(False)
