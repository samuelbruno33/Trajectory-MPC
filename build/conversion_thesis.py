import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ==== CONFIGURAZIONI ====
csv_path = "comparison_results.csv"  # Il file generato dal C++
plt.style.use('seaborn-v0_8-whitegrid') # Stile simile a quello di Fede

# ==== LETTURA FILE ====
df = pd.read_csv(csv_path)

# ==== 1. TRAIETTORIA (X vs Y) - Come Fig. 1 di Fede ====
# Confrontiamo MPC (che sterza) con PID (che va dritto/sbaglia)
plt.figure(figsize=(10, 6))
plt.plot(df['MPC_X'], df['MPC_Y'], label='Traiettoria MPC (Coupled)', linewidth=2, color='blue')
plt.plot(df['PID_X'], df['PID_Y'], label='Traiettoria PID (Baseline)', linestyle='--', color='red')
# Creiamo una finta traiettoria di riferimento per il plot (la sinusoide usata nel C++)
ref_x = df['MPC_X']
ref_y = 2.0 * np.sin(ref_x / 10.0)
plt.plot(ref_x, ref_y, label='Traiettoria di Riferimento', linestyle=':', color='yellow', linewidth=2)

plt.xlabel('X [m]', fontsize=12)
plt.ylabel('Y [m]', fontsize=12)
plt.title('Traiettoria: MPC vs PID Baseline', fontsize=14)
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.tight_layout()
plt.savefig('grafico_1_traiettoria.png', dpi=300)
plt.show()

# ==== 2. LONGITUDINALE: VELOCITÀ & POWER (Il Tuo contributo) ====
fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

# Velocità
axs[0].plot(df['Time'], df['Ref_Vx'], 'k--', label='Riferimento', linewidth=2)
axs[0].plot(df['Time'], df['MPC_Vx'], 'b-', label='MPC Velocity')
axs[0].plot(df['Time'], df['PID_Vx'], 'r-.', label='PID Velocity')
axs[0].set_ylabel('Velocità [m/s]', fontsize=12)
axs[0].set_title('Tracking Velocità Longitudinale', fontsize=14)
axs[0].legend()
axs[0].grid(True)

# Potenza (Constraint)
power_limit = 80000 # 80 kW
axs[1].plot(df['Time'], df['MPC_Fx'] * df['MPC_Vx'], 'b-', label='MPC Power')
axs[1].plot(df['Time'], df['PID_Power'], 'r-.', label='PID Power')
axs[1].axhline(y=power_limit, color='green', linestyle='-', linewidth=2, label='Limite 80kW')
# Evidenzia violazione PID
axs[1].fill_between(df['Time'], df['PID_Power'], power_limit, 
                    where=(df['PID_Power'] > power_limit), color='red', alpha=0.3, label='Violazione PID')

axs[1].set_xlabel('Tempo [s]', fontsize=12)
axs[1].set_ylabel('Potenza [W]', fontsize=12)
axs[1].set_title('Gestione Vincoli Potenza (MPC vs PID)', fontsize=14)
axs[1].legend()
axs[1].grid(True)

plt.tight_layout()
plt.savefig('grafico_2_longitudinale.png', dpi=300)
plt.show()

# ==== 3. STERZO E INPUT (Come Fig. 3 di Fede) ====
plt.figure(figsize=(10, 5))
plt.plot(df['Time'], df['MPC_Steer'], label='Sterzo MPC [rad]', color='blue')
# Il PID non ha sterzo nel nostro test (o è 0), lo mostriamo per confronto
plt.plot(df['Time'], np.zeros_like(df['Time']), label='Sterzo PID (Nullo)', color='red', linestyle='--')

plt.xlabel('Tempo [s]')
plt.ylabel('Angolo Sterzo [rad]')
plt.title('Attività di Sterzo (Controllo Laterale)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('grafico_3_sterzo.png', dpi=300)
plt.show()