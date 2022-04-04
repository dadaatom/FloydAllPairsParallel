#include <math.h>
#include <iostream>
#include <fstream>
#include <bitset>
#include <string>
#include <random>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"
#include <stdio.h>

int main(int argc, char* argv[]) {
	const int INF = 1000000;
	const int n = 36;

	std::vector<std::vector<int>> A = std::vector<std::vector<int>>();
	A.resize(36, std::vector<int>(36));

	int TempA[] = {
	  0,   1, INF, INF,   4, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF,   0,   2, INF,   3, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF,   0,   3,   2, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF,   0,   4, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF,   0,   5,   2, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,  20, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	  6, INF, INF, INF, INF,   0, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF,   0,   1, INF, INF, INF, INF,   6, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF,   0,   2, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF,   0,   3,   2, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF,   0,   4, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0,   5, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF,   6, INF,   3, INF, INF,   0,   1, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0,   1, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0,   2, INF,   3, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0,   3, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0,   4, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0,   5, INF, INF, INF, INF, INF,  10, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   6, INF, INF, INF, INF,   0, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0, INF,   2, INF, INF,   6, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   5,   0, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   4,   0, INF,   2, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   3,   0, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,  14,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   2,   0, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   2,   1,   0, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0,   1, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	 25, INF, INF,  10, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0,   2, INF, INF, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0,   1,   2, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0,   3, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   4,   3, INF, INF,   0, INF, INF, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   0, INF, INF, INF, INF,   6, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   2,   5,   0, INF, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   4,   0, INF, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   2,   3,   0, INF, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   2,   0, INF, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   3, INF,   1,   0, INF,
	INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,   1,   0 };

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = TempA[i*n+j];
		}
	}

	int id = 0;
	int p = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if (((int)std::sqrt(p)) == std::sqrt(p)) {
		int blockSize = n / std::sqrt(p);
		int row = id / std::sqrt(p) + 1;
		int col = id % (int)std::sqrt(p) + 1;

		MPI_Comm rowComm;
		MPI_Comm colComm;
		MPI_Comm_split(MPI_COMM_WORLD, row, id, &rowComm);
		MPI_Comm_split(MPI_COMM_WORLD, col, id, &colComm);

		/*
		int rowId;
		int colId;

		MPI_Comm_rank(rowComm, &rowId);
		MPI_Comm_rank(colComm, &colId);
		*/

		std::vector<std::vector<int>> G = std::vector<std::vector<int>>();
		G.resize(blockSize, std::vector<int>(blockSize));

		double parallelStart;
		if (id == 0) {
			parallelStart = MPI_Wtime();
			std::cout << "Initialising local matrices...\n";
		}

		for (int i = 0; i < blockSize; i++) {
			for (int j = 0; j < blockSize; j++) {
				int matrixI = (row - 1) * blockSize + i;
				int matrixJ = (col - 1) * blockSize + j;

				G[i][j] = A[matrixI][matrixJ];
			}
		}

		if (id == 0) {
			std::cout << "Executing iterations...\n";
		}

		int rootRow = 0;
		int rootCol = 0;

		for (int k = 0; k < n; k++) {
			if (k && ((k % blockSize) == 0)) {
				rootRow++;
				rootCol++;
			}

			//int rootRow = k / blockSize;
			//int rootCol = k / blockSize;

			//int rootId = rootRow * n + rootCol;

			std::vector<int> R = std::vector<int>(blockSize);
			std::vector<int> C = std::vector<int>(blockSize);

			if ((row - 1) * blockSize <= k && (row)*blockSize > k) {//rowId == ((int)(k / blockSize))) {//
				for (int i = 0; i < blockSize; i++) {
					R[i] = G[(k - (row - 1) * blockSize)][i];
				}
			}

			//BCAST ROW
			//MPI_Barrier(colComm);

			MPI_Bcast(&R[0], blockSize, MPI_INT, rootRow, colComm);

			if ((col - 1) * blockSize <= k && (col)*blockSize > k) {//colId == ((int)(k/blockSize))){//
				for (int i = 0; i < blockSize; i++) {
					C[i] = G[i][(k - (col - 1) * blockSize)];
				}
			}

			//BCAST COL
			MPI_Bcast(&C[0], blockSize, MPI_INT, rootCol, rowComm);

			//WAIT TO RECIEVE DATA.
			MPI_Barrier(rowComm);

			//COMPUTE PART OF THE MATRIX USING COLLECTED ROW/COL
			for (int i = 0; i < blockSize; i++) {
				for (int j = 0; j < blockSize; j++) {
					G[i][j] = std::min(G[i][j], C[i] + R[j]);
				}
			}
		}

		if (id == 0) {

			std::cout << "Collecting final matrices...\n";

			std::vector<std::vector<int>> H = std::vector<std::vector<int>>();
			H.resize(n, std::vector<int>(n));

			for (int i = 0; i < blockSize; i++) {
				for (int j = 0; j < blockSize; j++) {
					int matrixI = (row - 1) * blockSize + i;
					int matrixJ = (col - 1) * blockSize + j;

					H[matrixI][matrixJ] = G[i][j];
				}
			}

			for (int x = 1; x < p; x++) {
				MPI_Status status;

				std::vector<int> recv = std::vector<int>(blockSize * blockSize);
				MPI_Recv(&recv[0], blockSize * blockSize, MPI_INT, x, 0, MPI_COMM_WORLD, &status);

				int matrixRow = x / std::sqrt(p);
				int matrixCol = x % (int)std::sqrt(p);

				for (int i = 0; i < blockSize; i++) {
					for (int j = 0; j < blockSize; j++) {
						int matrixI = (matrixRow) * blockSize + i;
						int matrixJ = (matrixCol) * blockSize + j;

						H[matrixI][matrixJ] = recv[i * blockSize + j];
					}
				}
			}

			double parallelEnd = MPI_Wtime();

			std::cout << "Printing resulting parallel matrix...\n";
			for (int i = 0; i < n; i++) {
				std::cout << "\t[ ";
				for (int j = 0; j < n; j++) {
					std::cout << H[i][j] << " ";
				}
				std::cout << "]\n";
			}

			int serialStart = MPI_Wtime();

			for (int k = 0; k < n; k++) {
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n; j++) {
						if ((A[i][k] < INF && A[k][j] < INF)/* && (A[i][j] > A[i][k] + A[k][j])*/) {
							//A[i][j] = A[i][k] + A[k][j]; 
							A[i][j] = std::min(A[i][j], A[i][k] + A[k][j]);
						}
						/*if ((TempA[i*n + k] < INF && TempA[k * n + j] < INF) && (TempA[i * n + j] > TempA[i * n + k] + TempA[k * n + j])) {
							TempA[i * n + j] = TempA[i * n + k] + TempA[k * n + j]; //std::min(A[i][j], A[i][k] + A[k][j]);
						}*/
					}
				}
			}

			int serialEnd = MPI_Wtime();

			std::cout << "Printing resulting sequential matrix...\n";
			for (int i = 0; i < n; i++) {
				std::cout << "\t[ ";
				for (int j = 0; j < n; j++) {
					std::cout << A[i][j] << " ";
					//std::cout << TempA[i * n + j] << " ";
				}
				std::cout << "]\n";
			}
			
			std::cout << "Parallel process completed in " << (parallelEnd - parallelStart) << " seconds!\n";
			std::cout << "Serial process completed in " << (serialEnd - serialStart) << " seconds!\n";
			std::cout << "Speedup: " << ((serialEnd - serialStart) / (parallelEnd - parallelStart)) << "!\n";
			
		}
		else {
			std::vector<int> send = std::vector<int>(blockSize*blockSize);

			for (int i = 0; i < blockSize; i++) {
				for (int j = 0; j < blockSize; j++) {
					send[i * blockSize + j] = G[i][j];
				}
			}

			MPI_Send(&send[0], blockSize*blockSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}
	else {
		std::cout << "Number of processes is not a perfect square...\n";
	}

	MPI_Finalize();
	return 0;
}