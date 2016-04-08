Projet de Multicore_Programming
Auteurs : Léo Cassiau, Ugo Mahey

COMPILATION

Pour compiler le programme, éxecuter la commande ./make dans un terminal.
Trois éxecutables sont alors générés :
	- optimization-seq, qui le programme original qui cherche le minimum d'une fonction séquentiellement
	- optimization-omp, qui est optimization-seq modifié à l'aide d'OpenMP afin de paralléliser le programme sur une seule machine
	- optimization-mpi, qui est

EXECUTION

Pour éxecuter optimization-mpi : 
hostname 																					// Pour connaître le nom du pc.
mpirun -H I125V2pc14 -n 2 optimization-mpi 				// -H : indique l'hôte. -n : indique le nombre de machines.
