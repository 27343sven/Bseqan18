from SeqModules.Aligner import Aligner

if __name__ == '__main__':
    align = Aligner("actcgattgcct", "acttccgaatttggct")
    align.smith_waterman()