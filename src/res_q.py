


toFileName = "/Users/greg/Desktop/Residue_Qscores"

r = 2.9
mean=-0.0016*r*r*r+0.0434*r*r-0.3956*r+1.3366
peak=mean+0.024
low=mean-0.126
high=mean+0.109

print("Q - peak:%.2f, low:%.2f, high:%.2f"%(peak,low,high))

sel_residues = session.selection.items('residues')[0]
chain = sel_residues[0].chain_id


outf = open ( toFileName + "_%s.txt" % chain, "w" )
outf.write ( "Residue\tNum\tQ-score\tPeak\tLow_Q\tHigh_Q\n" )

for res in sel_residues :
    qscores = [at.qscore for at in res.atoms]
    res_q = sum (qscores) / len (qscores)
    outf.write ( "%s\t%d\t%f\t%f\t%f\t%f\n" % (res.name, res.number, res_q, peak, low, high ) )

outf.close()
