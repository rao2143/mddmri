function atomind = fMolnrAtomnam2Atomind(gro,molnr,atomnam)

atomind = find(all([gro.molnr==molnr strcmp(gro.atom,atomnam)],2));