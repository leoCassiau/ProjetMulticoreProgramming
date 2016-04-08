Projet de Multicore_Programming
Auteurs : Léo Cassiau, Ugo Mahey


COMPILATION

Pour compiler le programme, éxecuter la commande ./make dans un terminal.
Trois éxecutables sont alors générés :
	- optimization-seq, qui le programme original qui cherche le minimum d'une fonction séquentiellement
	- optimization-omp, qui est optimization-seq modifié à l'aide d'OpenMP afin de paralléliser le programme sur plusieurs processeurs d'une seule machine.
	- optimization-mpi, qui est optimization-omp modifié à l'aide de MPI afin de paralléliser le programme sur plusieurs processeurs de plusieurs machines.


EXECUTION

optimization-seq et optimization-omp s'éxecutent avec les commandes :
./optimization-seq
./optimization-omp

Pour éxecuter optimization-mpi : 
hostname	// Commande utile pour connaître le nom du pc.
mpirun -H I125V2pc13,I125V2pc14 -n 2 optimization-mpi	// Execute optimization-mpi sur les différentes machines indiqués. -H : indique les noms des machines hôtes. -n : indique le nombre de machines. Ici, I125V2pc13 et I125V2pc14 sont les 2 hôtes, I125V2pc13 étant la machine de rang 0 et I125V2pc14 de rang 1.

A la faculté des sciences et techniques de nantes, nous avons utilisé le mpirun de l'enseignant M. Goualard, ce qui donne la commande suivante :
~goualard-f/local/bin/mpirun -H I125V2pc13,I125V2pc14 -n 2 optimization-mpi
