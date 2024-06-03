function atomind = fMolnamAtomnam2Atomind(gro,molnam,atomnam)

atomind = find(all([strcmp(gro.molecule,molnam) strcmp(gro.atom,atomnam)],2));
