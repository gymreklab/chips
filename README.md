# asimon-chip-sim
Simulation tool for ChIP- and other -seq experiments


### Usage
asimon reads [options]

Required arguments:
1. -p <peaks>:
    The path to your peak file.
    Currently, ChIPMunk support BED format, and homer format. You need to specify the format with "-t".
2. -t <str>:
    The file format of your peak file.
    You can choose from "bed" and "homer".
