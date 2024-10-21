from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

a = SeqRecord(Seq("AAAACGT"), id="Alpha")
b = SeqRecord(Seq("AAA-CGT"), id="Beta")
c = SeqRecord(Seq("AAAAGGT"), id="Gamma")
align = MultipleSeqAlignment([a, b, c], annotations={"tool": "demo"}, column_annotations={"stats": "CCCXCCC"})

calculator = DistanceCalculator('identity')
distMatrix = calculator.get_distance(align)

constructor = DistanceTreeConstructor()

UPGMATree = constructor.upgma(distMatrix)
Phylo.draw(UPGMATree)
Phylo.draw_ascii(UPGMATree)
