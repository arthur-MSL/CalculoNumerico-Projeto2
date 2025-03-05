#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 9
#define EPSILON 1e-8

void salvar_csv(const char *metodo, int iteracao, double erro) {
    FILE *arquivo = fopen("resultados.csv", "a");
    if (!arquivo) {
        perror("Erro ao abrir arquivo");
        exit(EXIT_FAILURE);
    }
    fprintf(arquivo, "%s,%d,%.10f\n", metodo, iteracao, erro);
    fclose(arquivo);
}

void metodo_jacobi(double A[N][N], double B[N]) {
    double X[N] = {0};
    double X_novo[N];
    int iteracao = 0;
    double erro;
    
    do {
        erro = 0.0;
        for (int i = 0; i < N; i++) {
            double soma = 0.0;
            for (int j = 0; j < N; j++) {
                if (i != j) soma += A[i][j] * X[j];
            }
            X_novo[i] = (B[i] - soma) / A[i][i];
            erro = fmax(erro, fabs(X_novo[i] - X[i]));
        }
        for (int i = 0; i < N; i++) X[i] = X_novo[i];
        iteracao++;
        salvar_csv("Jacobi", iteracao, erro);
    } while (erro > EPSILON);
}

void metodo_gauss_seidel(double A[N][N], double B[N]) {
    double X[N] = {0};
    int iteracao = 0;
    double erro;
    
    do {
        erro = 0.0;
        for (int i = 0; i < N; i++) {
            double soma = 0.0;
            for (int j = 0; j < i; j++) soma += A[i][j] * X[j];
            for (int j = i + 1; j < N; j++) soma += A[i][j] * X[j];
            double novo_Xi = (B[i] - soma) / A[i][i];
            erro = fmax(erro, fabs(novo_Xi - X[i]));
            X[i] = novo_Xi;
        }
        iteracao++;
        salvar_csv("Gauss-Seidel", iteracao, erro);
    } while (erro > EPSILON);
}

void metodo_sor(double A[N][N], double B[N], double w, const char *nome) {
    double X[N] = {0};
    int iteracao = 0;
    double erro;
    
    do {
        erro = 0.0;
        for (int i = 0; i < N; i++) {
            double soma = 0.0;
            for (int j = 0; j < i; j++) soma += A[i][j] * X[j];
            for (int j = i + 1; j < N; j++) soma += A[i][j] * X[j];
            double novo_Xi = (1 - w) * X[i] + (w * (B[i] - soma) / A[i][i]);
            erro = fmax(erro, fabs(novo_Xi - X[i]));
            X[i] = novo_Xi;
        }
        iteracao++;
        salvar_csv(nome, iteracao, erro);
    } while (erro > EPSILON);
}

int main() {
    double A[N][N] = {
        {-4, 1, 0, 1, 0, 0, 0, 0, 0},
        {1, -4, 1, 0, 1, 0, 0, 0, 0},
        {0, 1, -4, 0, 0, 1, 0, 0, 0},
        {1, 0, 0, -4, 1, 0, 1, 0, 0},
        {0, 1, 0, 1, -4, 1, 0, 1, 0},
        {0, 0, 1, 0, 1, -4, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, -4, 1, 0},
        {0, 0, 0, 0, 1, 0, 1, -4, 1},
        {0, 0, 0, 0, 0, 1, 0, 1, -4}
    };
    double B[N] = {-50, -50, -150, 0, 0, -100, -50, -50, -150};

    FILE *arquivo = fopen("resultados.csv", "w");
    fprintf(arquivo, "Metodo,Iteracao,Erro_Final\n");
    fclose(arquivo);
    
    metodo_jacobi(A, B);
    metodo_gauss_seidel(A, B);
    metodo_sor(A, B, 0.5, "SOR_w=0.5");
    metodo_sor(A, B, 0.9, "SOR_w=0.9");
    metodo_sor(A, B, 1.2, "SOR_w=1.2");
    metodo_sor(A, B, 1.5, "SOR_w=1.5");
    
    printf("Resultados salvos em resultados.csv\n");
    return 0;
}
