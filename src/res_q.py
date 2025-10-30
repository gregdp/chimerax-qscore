

# change to local file path to output residue q-scores:
toFile = "/Users/greg/Desktop/Residue_Qscores"
# change to resolution of map:
r = 1.9

# !! change "Reference Sigma" in Q-score dialog to 0.4 !!

mean = -0.0016*r*r*r+0.0434*r*r-0.3956*r+1.3366
peak = mean+0.024
low = mean-0.126
high = mean+0.109
print("\nQ-score at %.2f A - peak:%.2f, low:%.2f, high:%.2f"%(r, peak,low,high))

sel_residues_items = session.selection.items('residues')
if len( sel_residues_items ) > 0 :
    sel_residues = sel_residues_items[0]
    if len(sel_residues) == 1 :
        res = sel_residues[0]
        qscores = [at.qscore for at in res.atoms if at.element != "H"]
        res_q = sum (qscores) / len (qscores)
        print ( "Residue: %s %s/:%d Q-score: %.2f\n" % (res.name, res.chain_id, res.number, res_q) )
    else :
        chain = sel_residues[0].chain_id
        outf = open ( toFile + "_%s.txt" % chain, "w" )
        outf.write ( "Residue\tNum\tQ-score\tPeak\tLow_Q\tHigh_Q\n" )
        for res in sel_residues :
            qscores = [at.qscore for at in res.atoms if at.element != "H"]
            res_q = sum (qscores) / len (qscores)
            outf.write ( "%s\t%d\t%f\t%f\t%f\t%f\n" % (res.name, res.number, res_q, peak, low, high ) )
        outf.close()
