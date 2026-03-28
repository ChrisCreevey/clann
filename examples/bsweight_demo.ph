[
 Bootstrap-weighted supertree demonstration dataset
 ===================================================
 8 taxa: Human Chimp Gorilla Orang Mouse Rat Dog Cat

 TRUE topology:      (((Human,Chimp),(Gorilla,Orang)),((Mouse,Rat),(Dog,Cat)))
 COMPETING topology: (((Human,Gorilla),(Chimp,Orang)),((Mouse,Rat),(Dog,Cat)))

 The two topologies differ only in the primate groupings:
   TRUE:      (Human+Chimp) and (Gorilla+Orang) as separate sister pairs
   COMPETING: (Human+Gorilla) and (Chimp+Orang) as sister pairs

 The dataset contains:
   4 trees supporting the TRUE topology      — HIGH bootstrap support (92-99)
                                               on the primate grouping nodes
   7 trees supporting the COMPETING topology — LOW  bootstrap support (27-41)
                                               on the same primate grouping nodes

 Expected behaviour
 ------------------
 WITHOUT bsweight (default): the 7-vs-4 majority vote favours the COMPETING
   topology, which is incorrect relative to the high-confidence signal.

 WITH bsweight=1: the low-confidence primate splits in the competing trees are
   downweighted proportionally to their bootstrap support. The 4 high-confidence
   trees supporting the TRUE topology outweigh the 7 uncertain trees, recovering
   the correct grouping.

 Suggested commands
 ------------------
   exe bsweight_demo.ph
   set criterion sfit

   hs nreps=10 progress=1           (bsweight off — majority wins)
   hs nreps=10 progress=1 bsweight 1   (bsweight on  — confidence wins)
]

[TRUE topology, HIGH bootstrap support — 4 trees]
(((Human,Chimp)95,(Gorilla,Orang)96)98,((Mouse,Rat)97,(Dog,Cat)95)99);
(((Human,Chimp)93,(Gorilla,Orang)94)97,((Mouse,Rat)96,(Dog,Cat)93)98);
(((Human,Chimp)96,(Gorilla,Orang)95)99,((Mouse,Rat)98,(Dog,Cat)96)97);
(((Human,Chimp)92,(Gorilla,Orang)91)95,((Mouse,Rat)95,(Dog,Cat)94)97);

[COMPETING topology, LOW bootstrap support on primate nodes — 7 trees]
(((Human,Gorilla)32,(Chimp,Orang)35)38,((Mouse,Rat)97,(Dog,Cat)95)99);
(((Human,Gorilla)28,(Chimp,Orang)31)33,((Mouse,Rat)96,(Dog,Cat)93)98);
(((Human,Gorilla)35,(Chimp,Orang)38)41,((Mouse,Rat)98,(Dog,Cat)96)97);
(((Human,Gorilla)30,(Chimp,Orang)33)36,((Mouse,Rat)97,(Dog,Cat)94)99);
(((Human,Gorilla)27,(Chimp,Orang)29)32,((Mouse,Rat)95,(Dog,Cat)95)98);
(((Human,Gorilla)33,(Chimp,Orang)36)39,((Mouse,Rat)96,(Dog,Cat)92)97);
(((Human,Gorilla)31,(Chimp,Orang)34)37,((Mouse,Rat)97,(Dog,Cat)93)98);
