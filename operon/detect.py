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
# Roback, P et al. “A predicted operon map for Mycobacterium tuberculosis.”
# Nucleic acids research vol. 35,15 (2007): 5085-95. doi:10.1093/nar/gkm518
from __future__ import division
from dataclasses import dataclass, field
from . import Analyze, average_coverage
from .gtf_process import GtfProcess
from .operon import Operon
import os.path
import pysam
import csv

@dataclass
class Detect(Analyze):

    ref_seq: str
    ref_len: int
    bam_file: str
    gtf_file: str
    gene_depth: int = field(default=2)
    igr_depth: int = field(default=1)
    gene_factor: float = field(default=5)
    igr_factor: float = field(default=5)
    output: str = field(default='detected-operons.csv')

    @classmethod
    def print(cls, fs_writer, op):
        """
        | | Write op info and op genes in csv file
        """
        row = op.info_csv_row()
        fs_writer.writerow(row)
        rows = op.gene_csv_rows()
        fs_writer.writerows(rows)

    def adjacent_coverages(self, op, avg_cov):
        """
        | | Compare operon gene expression coverages to the current gene coverage
        """
        adj_covs = []
        max_cov = max(op.genes_covs())
        min_cov = min(op.genes_covs())
        if avg_cov >= max_cov:
            val = min_cov
        elif avg_cov > min_cov:
            if max_cov / avg_cov > avg_cov / min_cov:
                val = max_cov
            else:
                val = min_cov
        else:
            val = max_cov


        for i, cov in enumerate(op.genes_covs()):
            #if cov == 0 or avg_cov ==0: #(Tracey added this line and the next and changes the third line after this to "elif". check Hocine)
                 #continue
            try:
            	if avg_cov / cov < self.gene_factor and cov / avg_cov < self.gene_factor:
                	adj_covs.append(cov)
            	else:
                	adj_covs.append(cov)
                	if cov == val:
                	    break
            except ZeroDivisionError:
                pass

        return adj_covs

    def iql_coverage(self, alignment_file, start, op_end, strand, cyclic=False):
        """
        | | Calculate the coverage of the interquartile region of IGR
        """
        if cyclic:
            mean = (start + (self.ref_len - op_end)) / 2
            first_iql = int(round(op_end + mean /2 ))
            last_iql = int(round(start + self.ref_len - mean / 2))
            if first_iql > self.ref_len:
                first_iql = first_iql - self.ref_len
                last_iql = start - mean / 2
                iql_cov = average_coverage(alignment_file, self.ref_seq, first_iql, last_iql, strand)
            elif last_iql > self.ref_len:
                cov_1 = average_coverage(alignment_file, self.ref_seq, first_iql, 0, strand)
                cov_2 = average_coverage(alignment_file, self.ref_seq, 0, last_iql - self.ref_len, strand)
                iql_cov = (cov_1 + cov_2) / 2
            else:
                iql_cov = average_coverage(alignment_file, self.ref_seq, first_iql, last_iql, strand)
        else:
            mean = (start - op_end) / 2
            first_iql = int(round(op_end + mean / 2))
            last_iql = int(round(start - mean / 2))
            iql_cov = average_coverage(alignment_file, self.ref_seq, first_iql, last_iql, strand)

        return iql_cov

    def analyze(self):
        """
        | |
        """
        # Create the output directory if it does not already exist.
        os.makedirs('./output', exist_ok=True)
        fs = open(os.path.join('./output/', self.output), 'w')
        fs_writer = csv.writer(fs, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        gtf_process = GtfProcess(self.gtf_file)
        alignment_file = pysam.AlignmentFile(self.bam_file, 'rb')#, index_filename='test/data/test_file.bai')
        lines = gtf_process.lines()

        op = None
        num_genes = 0
        for line in lines:
            start = int(line['start'])
            end = int(line['end'])
            avg_cov = average_coverage(alignment_file, self.ref_seq, start, end, line['strand'])
            if op is None:
                if avg_cov < self.gene_depth:
                    fs_writer.writerow([
                        line['gene_id'],
                        str(start),
                        str(end),
                        str(end - start),
                        line['strand'],
                        "NOT EXPRESSED",
                        str(avg_cov)
                    ])
                else:
                    op = Operon.operon(start, end, line['strand'], line['gene_id'], avg_cov)
            elif avg_cov < self.gene_depth:
                self.print(fs_writer, op)
                fs_writer.writerow([
                    line['gene_id'],
                    str(start),
                    str(end),
                    str(end - start),
                    line['strand'],
                    "NOT EXPRESSED",
                    str(avg_cov)
                ])
                op = None
            elif op.strand != line['strand']:
                self.print(fs_writer, op)
                op = Operon.operon(start, end, line['strand'], line['gene_id'], avg_cov)
            else:
                if op.end < start and start - op.end > 4:
                    igr_cov = average_coverage(alignment_file, self.ref_seq, op.end, start, line['strand'])
                    # mean = (start - op.end) / 2
                    # first_iql = int(round(op.end + mean / 2))
                    # last_iql = int(round(start - mean / 2))
                    # iql_cov = average_coverage(alignment_file, self.ref_seq, first_iql, last_iql, line['strand'])
                    iql_cov = self.iql_coverage(alignment_file, start, op.end, line['strand'])
                    last_gene_cov = op.last_gene_cov()
                    if igr_cov > self.igr_depth and iql_cov > self.igr_depth and \
                       last_gene_cov / iql_cov < self.igr_factor and iql_cov / last_gene_cov < self.igr_factor and \
                       avg_cov / iql_cov < self.igr_factor and iql_cov / avg_cov < self.igr_factor:

                        adj_covs = self.adjacent_coverages(op, avg_cov)
                        if len(adj_covs) == len(op.genes_covs()):
                            op.add_gene(start, end, line['gene_id'], avg_cov, igr_cov)
                        else:
                            if len(adj_covs) > 0:
                                op_1 = None
                                op_2 = None
                                is_diff = False
                                op_covs = op.genes_covs()
                                for i, gene in enumerate(op.genes):
                                    if i < len(adj_covs):
                                        if op_1 is None:
                                            op_1 = Operon.operon(op.start, op.end, op.strand, gene.gene_id, gene.cov)
                                        else:
                                            op_1.add_gene(gene.start, gene.end, gene.gene_id, gene.cov, gene.igr_cov)
                                    else:
                                        if not is_diff and op_1 is not None:
                                            r1 = adj_covs[-1] / op_covs[i] if op_covs[i] <= adj_covs[-1] else \
                                                 op_covs[i] / adj_covs[-1]
                                            r2 = avg_cov / op_covs[i] if op_covs[i] <= adj_covs[-1] else \
                                                 op_covs[i] / avg_cov
                                            if r1 <= r2:
                                                op_1.add_gene(gene.start, gene.end, gene.gene_id, gene.cov, gene.igr_cov)
                                                adj_covs.append(op_covs[i])
                                            else:
                                                self.print(fs_writer, op_1)
                                                op_2 = Operon.operon(gene.start, gene.end, op.strand, gene.gene_id, gene.cov)
                                                is_diff = True
                                        else:
                                            if is_diff:
                                                op_2.add_gene(gene.start, gene.end, gene.gene_id, gene.cov, gene.igr_cov)
                                                # op_2.genes.append(gene)
                                if op_2 is None and op_1 is not None:
                                    self.print(fs_writer, op_1)
                                    op_2 = Operon.operon(start, end, line['strand'], line['gene_id'], avg_cov)
                                elif op_2 is not None:
                                    gene = op_2.genes[-1]
                                    if gene.cov / avg_cov < self.gene_factor and avg_cov / gene.cov < self.gene_factor:
                                        op_2.add_gene(start, end, line['gene_id'], avg_cov, igr_cov)
                                    else:
                                        self.print(fs_writer, op_2)
                                        op_2 = Operon.operon(start, end, line['strand'], line['gene_id'], avg_cov)
                                op = op_2
                            else:
                                self.print(fs_writer, op)
                                op = Operon.operon(start, end, line['strand'], line['gene_id'], avg_cov)
                    else:
                        self.print(fs_writer, op)
                        op = Operon.operon(start, end, line['strand'], line['gene_id'], avg_cov)
                elif op.end < start:
                    igr_cov = average_coverage(alignment_file, self.ref_seq, op.end, start, line['strand'])
                    if igr_cov > self.igr_depth:
                        op.add_gene(start, end, line['gene_id'], avg_cov, igr_cov)
                else:
                    op.add_gene(start, end, line['gene_id'], avg_cov, 0)
            if num_genes % 1000 == 0:
                print('\n' + str(num_genes / 1000) + ':', end='')
            else:
                if num_genes % 100 == 0:
                    print('.', end=' ')
            num_genes += 1

        # Cyclic genome. continue with the first cds
        lines = gtf_process.lines()
        line_0 = lines.__next__()

        if op is not None:
            start = int(line_0['start'])
            end = int(line_0['end'])
            strand = line_0['strand']
            avg_cov = average_coverage(alignment_file, self.ref_seq, start, end, strand)
            if op.strand != strand or avg_cov < self.gene_depth:
                self.print(fs_writer, op)
                op = None
            else:
                # whole IGR coverage
                cov_1 = average_coverage(alignment_file, self.ref_seq, op.end, 0, strand)
                cov_2 = average_coverage(alignment_file, self.ref_seq, 0, start, strand)
                igr_cov = (cov_1 + cov_2) / 2
                if start + (self.ref_len - op.end) > 4:
                    iql_cov = self.iql_coverage(alignment_file, start, op.end, strand, True)
                    last_gene_cov = op.last_gene_cov()
                    if igr_cov > self.igr_depth and iql_cov > self.igr_depth and \
                       last_gene_cov / iql_cov < self.igr_factor and iql_cov / last_gene_cov < self.igr_factor and \
                       avg_cov / iql_cov < self.igr_factor and iql_cov / avg_cov < self.igr_factor:

                        op.add_gene(start, end, line_0['gene_id'], avg_cov, igr_cov)
                    else :
                        self.print(fs_writer, op)
                        op = None
                else:
                    if igr_cov > self.igr_depth:
                        op.add_gene(start, end, line_0['gene_id'], avg_cov, igr_cov)
                    else:
                        self.print(fs_writer, op)
                        op = None

        # Cyclic genome continue with the remaining cds
        num_genes = 1
        if op is not None:
            for line in lines:
                start = int(line['start'])
                end = int(line['end'])
                avg_cov = average_coverage(alignment_file, self.ref_seq, start, end, line['strand'])
                if avg_cov < self.gene_depth:
                    self.print(fs_writer, op)
                    break
                elif op.strand != line['strand']:
                    self.print(fs_writer, op)
                    break
                else:
                    igr_cov = average_coverage(alignment_file, self.ref_seq, start, end, line['strand'])
                    if op.end < start and (start - op.end) > 4:
                        last_gene_cov = op.last_gene_cov()
                        iql_cov = self.iql_coverage(alignment_file, start, op.end, line['strand'])
                        if igr_cov > self.igr_depth and iql_cov > self.igr_depth and \
                                last_gene_cov / iql_cov < self.igr_factor and iql_cov / last_gene_cov < self.igr_factor and \
                                avg_cov / iql_cov < self.igr_factor and iql_cov / avg_cov < self.igr_factor:
                            op.add_gene(start, end, line['gene_id'], avg_cov, igr_cov)
                        else:
                            self.print(fs_writer, op)
                            break
                    elif op.end < start:
                        if igr_cov > self.igr_depth:
                            op.add_gene(start, end, line['gene_id'], avg_cov, igr_cov)
                        else:
                            self.print(fs_writer, op)
                            break
                    else:
                        op.add_gene(start, end, line['gene_id'], avg_cov, 0)

                if num_genes % 1000 == 0:
                    print('\n' + str(num_genes / 1000) + ':', end='')
                else:
                    if num_genes % 100 == 0:
                        print('.', end=' ')
                num_genes += 1

        if op is not None:
            self.print(fs_writer, op)
            op = None

        alignment_file.close()
