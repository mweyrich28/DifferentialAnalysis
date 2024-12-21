import matplotlib.pyplot as plt
from typing import List, Tuple, Optional, Union, Dict
import numpy as np
from collections import defaultdict


class GenomicIntervalPlotter:
    def __init__(self, figure_size: Tuple[int, int] = (12, 6),
                 exon_color: str = 'red',
                 interval_color: str = 'blue',
                 fw_color: str = 'darkgreen',
                 rw_color: str = 'purple',
                 font_size: int = 10,
                 vertical_spacing: float = 0.2):
        self.figure_size = figure_size
        self.exon_color = exon_color
        self.interval_color = interval_color
        self.fw_color = fw_color
        self.rw_color = rw_color
        self.font_size = font_size
        self.vertical_spacing = vertical_spacing

    @staticmethod
    def _parse_interval(interval: str) -> Tuple[int, int]:
        """Parse interval string and validate the format."""
        try:
            start, stop = map(int, interval.split('-'))
            if start > stop:
                raise ValueError(f"Invalid interval: start ({start}) > stop ({stop})")
            return start, stop
        except ValueError as e:
            raise ValueError(f"Invalid interval format: {interval}. Expected 'start-stop'") from e

    @staticmethod
    def _parse_read_with_id(read: str) -> Tuple[str, Tuple[int, int], str]:
        """Parse read string with ID, interval, and direction."""
        try:
            parts = read.split()
            if len(parts) != 3:
                raise ValueError("Expected format: 'id: start-stop direction'")
            
            read_id_part, direction = parts[0], parts[2]
            read_id = read_id_part.rstrip(':')
            interval = parts[1]
            
            if direction not in ['FW', 'RW']:
                raise ValueError(f"Invalid direction: {direction}. Expected 'FW' or 'RW'")
            
            start, stop = map(int, interval.split('-'))
            if start > stop:
                raise ValueError(f"Invalid interval: start ({start}) > stop ({stop})")
                
            return read_id, (start, stop), direction
            
        except ValueError as e:
            raise ValueError(f"Invalid read format: {read}. Expected 'id: start-stop direction'") from e

    def plot_intervals(self, 
                      exons: List[str], 
                      gaps: List[str],
                      reads: List[str],
                      title: Optional[str] = "Genomic Intervals",
                      show_coordinates: bool = True,
                      min_padding: int = 100) -> None:

        # Create figure
        plt.figure(figsize=self.figure_size)
        
        # Group reads by ID
        read_groups = defaultdict(lambda: {'FW': [], 'RW': []})
        all_positions = []
        
        # Process exons and gaps
        for interval_list in [exons, gaps]:
            for interval in interval_list:
                start, stop = self._parse_interval(interval)
                all_positions.extend([start, stop])
        
        # Process reads and group them by ID
        for read in reads:
            read_id, (start, stop), direction = self._parse_read_with_id(read)
            read_groups[read_id][direction].append((start, stop))
            all_positions.extend([start, stop])
        
        min_pos, max_pos = min(all_positions), max(all_positions)
        padding = max(min_padding, (max_pos - min_pos) * 0.05)
        
        # Plot exons
        for idx, exon in enumerate(exons):
            start, stop = self._parse_interval(exon)
            y_pos = 0  # Fixed position for exons
            plt.plot([start, stop], [y_pos, y_pos], 
                    color=self.exon_color, 
                    linewidth=4, 
                    label=f"Exon {idx + 1}" if idx == 0 else "")
            
            if show_coordinates:
                plt.text(start - padding/10, y_pos + 0.05, f"{start:,}", 
                        color='black', fontsize=self.font_size, ha='right', rotation=90)
                plt.text(stop + padding/10, y_pos - 0.07, f"{stop:,}", 
                        color='black', fontsize=self.font_size, ha='left', rotation=90)

        # Plot gaps
        for idx, gap in enumerate(gaps):
            y_pos = (idx + 1) * self.vertical_spacing
            
            gap_start, gap_stop = self._parse_interval(gap)
            plt.plot([gap_start, gap_stop], [y_pos, y_pos], 
                    color=self.interval_color, 
                    linewidth=4, 
                    label=f"Gap {idx + 1}" if idx == 0 else "")
            
            if show_coordinates:
                plt.text(gap_start - padding/10, y_pos + 0.05, f"{gap_start:,}", 
                        color='black', fontsize=self.font_size, ha='right', rotation=90)
                plt.text(gap_stop + padding/10, y_pos - 0.07, f"{gap_stop:,}", 
                        color='black', fontsize=self.font_size, ha='left', rotation=90)

        # Plot read groups
        for idx, (read_id, directions) in enumerate(read_groups.items()):
            y_pos = (len(gaps) + idx + 1) * self.vertical_spacing
            
            # Combine and sort all intervals while keeping direction information
            all_intervals = [(start, stop, 'FW') for start, stop in sorted(directions['FW'])]
            all_intervals.extend([(start, stop, 'RW') for start, stop in sorted(directions['RW'])])
            all_intervals.sort(key=lambda x: x[0])  # Sort by start position
            
            # Plot intervals and connections
            for i, (start, stop, direction) in enumerate(all_intervals):
                # Plot interval
                color = self.fw_color if direction == 'FW' else self.rw_color
                plt.plot([start, stop], [y_pos, y_pos], 
                        color=color, 
                        linewidth=4, 
                        label=f"Read {read_id} {direction}" if i == 0 else "")
                
                # Add coordinates if needed
                if show_coordinates:
                    plt.text(start - padding/10, y_pos + 0.05, f"{start:,}", 
                            color='black', fontsize=self.font_size, ha='right', rotation=90)
                    plt.text(stop + padding/10, y_pos - 0.07, f"{stop:,}", 
                            color='black', fontsize=self.font_size, ha='left', rotation=90)
                
                # Add connecting line to next interval if it exists
                if i < len(all_intervals) - 1:
                    next_start = all_intervals[i + 1][0]
                    next_direction = all_intervals[i + 1][2]
                    
                    # Choose line style based on direction comparison
                    if direction == next_direction:
                        # Same direction: solid line
                        linestyle = '-'
                    else:
                        # Different direction: dotted line
                        linestyle = ':'
                    
                    plt.plot([stop, next_start], [y_pos, y_pos],
                            color=color,
                            linewidth=1,
                            linestyle=linestyle,
                            alpha=0.5)

        # Customize plot
        plt.xlim(min_pos - padding, max_pos + padding)
        plt.xlabel("Genomic Position", fontsize=self.font_size + 2)
        plt.ylabel("Tracks", fontsize=self.font_size + 2)
        plt.title(title, fontsize=self.font_size + 4, pad=20)
        plt.yticks([])
        
        # Add grid and legend
        plt.grid(axis='x', linestyle='--', alpha=0.7)
        plt.legend(loc='upper right', bbox_to_anchor=(1, 1.15))
        
        # Adjust layout and display
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    # Create plotter with custom settings
    plotter = GenomicIntervalPlotter(
           figure_size=(12, 6),
           exon_color='red',
           interval_color='navy',
           fw_color='darkgreen',
           rw_color='purple',
           font_size=10,
           vertical_spacing=0.4
       )

    #   // Sample 3, 4, 6 wrong incl count (21 vs 20) → also plotted this case...
    exons = ["30702432-30702470"]
    gaps = []
    reads = [
        "66592: 30703331-30703402 FW",
        "66592: 30702122-30702169 RW",
        "66592: 30702001-30702033 RW",
        "66592: 30703584-30703611 FW",
    ]

    plotter.plot_intervals(exons, gaps, reads)

    # // Sample 6 wrong incl count (56 vs 55) → also plotted this case...
    exons = ["557720-557758"]
    gaps = []
    reads = ["36778: 555924-556023 RW", "36778: 563427-563526 FW"]

    plotter.plot_intervals(exons, gaps, reads)

    # // Sample 7
    # exons = ["30702432-30702471"]
    exons = [
        "30695943-30696118",
        "30696204-30696359",
        "30696455-30696552",
        "30696635-30696735",
        "30697310-30697401",
        "30697833-30697908",
        "30699154-30699222",
        "30699478-30699734",
        "30700111-30700195",
        "30700391-30700498",
        "30700605-30700691",
        "30701744-30701853",
        "30701949-30702033",
        "30702122-30702169",
        "30702432-30702470",
        "30703322-30703402",
        "30703584-30703643",
        "30704532-30704651",
        "30706318-30706365"
            ]
    gaps = []
    reads = [
        "82151: 30701980-30702033 RW",
        "82151: 30702122-30702167 RW",
        "82151: 30703323-30703402 FW",
        "82151: 30703584-30703588 FW",
    ]

    # Generate plot
    plotter.plot_intervals(exons, gaps, reads)
