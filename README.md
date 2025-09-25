# Distributed Matrix Multiplication using MPI (MapReduce Style)

## 📌 Overview
This project implements **parallel matrix multiplication** using the **Message Passing Interface (MPI)**.  
It follows a **MapReduce-like architecture**:
- **Mapper**: Converts matrix entries into key-value pairs.
- **Shuffler**: Rearranges keys to align multiplication pairs.
- **Reducer**: Performs partial multiplication and accumulation.
- **Gather**: Collects results and forms the final output matrix.

The project consists of two C programs:
1. **writing.c** → Generates random matrices A and B (stored in `MatrixA.txt` and `MatrixB.txt`).
2. **project.c** → Performs distributed matrix multiplication using MPI.

---

## 📂 File Structure
├── writing.c # Generates random matrices

├── project.c # MPI-based matrix multiplication

├── MatrixA.txt # Auto-generated matrix A (input)

├── MatrixB.txt # Auto-generated matrix B (input)

├── MatrixC.txt # Final result matrix (output)

└── README.md # Documentation


---

## ⚙️ Requirements
- **C Compiler** (e.g., `gcc`)
- **MPI Implementation** (e.g., [OpenMPI](https://www.open-mpi.org/) or MPICH)
- Linux/macOS (recommended). Windows users can use **WSL2** or Cygwin.

---

## 🔧 Installation
### 1. Install OpenMPI (Ubuntu/Debian)
```bash
sudo apt update
sudo apt install -y build-essential openmpi-bin libopenmpi-dev
```

### 2. Verify Installation
mpicc --version


## ▶️ How to Run
### Step 1: Generate Random Matrices
mpicc writing.c -o writing
./writing


This creates:

MatrixA.txt

MatrixB.txt

### Step 2: Run Parallel Matrix Multiplication
mpicc project.c -o project -lm
mpirun -np 4 ./project


-np 4 means run with 4 processes (1 master + 3 workers).

Adjust number of processes based on CPU cores available.

### Step 3: Check Results

Output matrix is written to MatrixC.txt.

The program also checks correctness by comparing with sequential multiplication.

## 📜 License

This project is open-source and free to use for learning and research purposes.