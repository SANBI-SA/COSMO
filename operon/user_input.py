# Copyright (C) 2019  Hocine Bendou <hocine@sanbi.ac.za>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>
from .detect import Detect
import argparse


def run_command(args: argparse.ArgumentParser):
    """
    Execute user's input command and return detected operons in a csv file
    """
    detect = Detect(args.ref, args.length, args.bam, args.gtf, args.gdepth, args.idepth,
                    args.gfactor, args.ifactor, args.output)
    detect.analyze()

    
def main():

    parser = argparse.ArgumentParser(description="Detect possible genome operons using RNA expression coverages")

    # positional arguments
    parser.add_argument("ref", help="Name of reference sequence in BAM file", type=str)
    parser.add_argument("length", help="Length of reference sequence", type=int)
    parser.add_argument("bam", help="Bam input file", type=argparse.FileType())
    parser.add_argument("gtf", help="Gtf input file")

    # optional arguments
    parser.add_argument("-D", "--gdepth", help="Average number of reads per base required to consider a gene expressed",
                                          type=int, default=2)
    parser.add_argument("-d", "--idepth", help="Average number of reads per base required to consider a IGR expressed",
                                          type=int, default=1)
    parser.add_argument("-F", "--gfactor", help="Allowed difference factor of two gene's coverages to be part of "
                                                "same Operon", type=float, default=5.0)
    parser.add_argument("-f", "--ifactor", help="Allowed difference factor of IGR and adjacent gene's coverages "
                                                "to be part of same operon", type=float, default=5.0)
    parser.add_argument("-o", "--output", help="CSV output filename under output folder", type=str, default="detected-operons.csv")

    args = parser.parse_args()
    run_command(args)

if __name__ == "__main__":
    main()
