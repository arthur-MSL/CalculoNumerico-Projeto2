import pandas as pd
import matplotlib.pyplot as plt


caminho_csv = "C:\\Users\\Arthu\\OneDrive\\Documentos\\UFU\\CNumerico\\projeto2\\resultados.csv"


df = pd.read_csv(caminho_csv)


plt.figure(figsize=(10, 6))

for metodo in df["Metodo"].unique():
    subset = df[df["Metodo"] == metodo]
    plt.plot(subset["Iteracao"], subset["Erro_Final"], marker='o', linestyle='-', label=metodo)

plt.xlabel("Iterações")
plt.ylabel("Erro Absoluto (norma ∞)")
plt.yscale("log") 
plt.title("Convergência dos Métodos Iterativos")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.xticks(range(0, df["Iteracao"].max() + 1, max(1, df["Iteracao"].max() // 10)))  
plt.show()
