import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ==== CONFIGURAZIONI ====
csv_path = "comparison_results.csv"
plt.style.use('seaborn-v0_8-whitegrid')

# ==== LETTURA FILE ====
df = pd.read_csv(csv_path)

# ==== 1. TRAIETTORIA (CHICANE TEST) ====
plt.figure(figsize=(12, 5))

# Plot Riferimento (Linea tratteggiata nera o gialla)
plt.plot(df['Ref_X'], df['Ref_Y'], 'k--', label='Reference Path', linewidth=1.5, alpha=0.6)

# Plot Veicoli
plt.plot(df['MPC_X'], df['MPC_Y'], label='MPC (Coupled)', linewidth=2.5, color='#1f77b4') # Blu
plt.plot(df['PID_X'], df['PID_Y'], label='PID + PurePursuit', linestyle='-.', linewidth=2, color='#d62728') # Rosso

plt.xlabel('Posizione X [m]', fontsize=12)
plt.ylabel('Posizione Y [m]', fontsize=12)
plt.title('Traiettoria: Chicane Test (High Speed Entry)', fontsize=14)
plt.legend()
plt.grid(True)
plt.axis('equal') # Importante per non distorcere la curva
plt.tight_layout()
plt.savefig('grafico_1_traiettoria.png', dpi=300)
print("Generato grafico_1_traiettoria.png")

# ==== 2. LONGITUDINALE COMPLETO (VELOCITA + POTENZA) ====
fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

# A. Velocità
axs[0].plot(df['Time'], df['Ref_Vx'], 'k--', label='Target Speed', linewidth=2)
axs[0].plot(df['Time'], df['MPC_Vx'], 'b-', label='MPC Speed')
axs[0].plot(df['Time'], df['PID_Vx'], 'r-.', label='PID Speed')
axs[0].set_ylabel('Velocità [m/s]', fontsize=12)
axs[0].set_title('Tracking Velocità (Frenata Rigenerativa & Ripresa)', fontsize=14)
axs[0].legend(loc='upper right')
axs[0].grid(True)

# B. Potenza (Il punto chiave della tesi)
power_limit = 30000 # 30 kW per il test
axs[1].plot(df['Time'], df['MPC_Fx'] * df['MPC_Vx'], 'b-', label='MPC Power Consumed')
axs[1].plot(df['Time'], df['PID_Power'], 'r-.', label='PID Power Req.')
axs[1].axhline(y=power_limit, color='green', linestyle='-', linewidth=2, label='Power Limit (30kW)')

# Evidenzia violazione PID
axs[1].fill_between(df['Time'], df['PID_Power'], power_limit, 
                    where=(df['PID_Power'] > power_limit), color='red', alpha=0.3, label='Violazione PID')

axs[1].set_xlabel('Tempo [s]', fontsize=12)
axs[1].set_ylabel('Potenza [W]', fontsize=12)
axs[1].set_title('Gestione Vincoli Potenza (Safety Constraint)', fontsize=14)
axs[1].legend(loc='upper right')
axs[1].grid(True)

plt.tight_layout()
plt.savefig('grafico_2_longitudinale.png', dpi=300)
print("Generato grafico_2_longitudinale.png")

plt.show()