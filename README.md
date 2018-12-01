# frame-9937-transKE-bug-pre-v1.0.0

bug in the calculation of trans_ke between frame 9937 and 9938. Bug appears only when using getFrameRangeEnergy() 
function but it does not applear while calculating separatly.
Possible cause is repeated use of shiftOrigin() twice sequencially; first in getRingRotKE() and second in getRingTransKE().
This bug could be eliminated by avoiding second use of shiftOrigin() in getRingTransKE()
