import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

csv_path = "comparison_results.csv"
plt.style.use('seaborn-v0_8-whitegrid')

# Lettura dati
df = pd.read_csv(csv_path)

# ============================================================
# GRAFICO 1: TRAIETTORIA
# ============================================================
plt.figure(figsize=(14, 6))

plt.plot(df['Ref_X'], df['Ref_Y'], 'k--', label='Traiettoria Riferimento', linewidth=1.5, alpha=0.6)
plt.plot(df['MPC_X'], df['MPC_Y'], label='MPC Accoppiato', linewidth=2.5, color='#1f77b4')
plt.plot(df['PID_X'], df['PID_Y'], label='PID + PurePursuit', linestyle='-.', linewidth=2, color='#d62728')

plt.xlabel('Posizione X [m]', fontsize=13)
plt.ylabel('Posizione Y [m]', fontsize=13)
plt.title('Confronto Traiettoria: MPC vs PID+PurePursuit', fontsize=15, fontweight='bold')
plt.legend(loc='upper right', fontsize=11)
plt.grid(True, alpha=0.3)
plt.axis('equal')
plt.tight_layout()
plt.savefig('grafico_1_traiettoria.png', dpi=300, bbox_inches='tight')
print("Generato: grafico_1_traiettoria.png")

# ============================================================
# GRAFICO 2: VELOCITA' E POTENZA
# ============================================================
fig, axs = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

# Subplot velocità
axs[0].plot(df['Time'], df['Ref_Vx'], 'k--', label='Riferimento', linewidth=2.5)
axs[0].plot(df['Time'], df['MPC_Vx'], 'b-', label='MPC', linewidth=2)
axs[0].plot(df['Time'], df['PID_Vx'], 'r-.', label='PID', linewidth=2)
axs[0].set_ylabel('Velocità [m/s]', fontsize=13)
axs[0].set_title('Tracking Velocità Longitudinale', fontsize=15, fontweight='bold')
axs[0].legend(loc='upper right', fontsize=11, framealpha=0.9)
axs[0].grid(True, alpha=0.3)
axs[0].set_ylim([0, 28])

# Subplot potenza (solo valori positivi)
mpc_power_motor = np.maximum(0, df['MPC_Fx'] * df['MPC_Vx'])
pid_power_motor = np.maximum(0, df['PID_Fx'] * df['PID_Vx'])
power_limit = 30000

axs[1].plot(df['Time'], mpc_power_motor, 'b-', label='MPC Potenza Motore', linewidth=2)
axs[1].plot(df['Time'], pid_power_motor, 'r-.', label='PID Potenza Motore', linewidth=2)
axs[1].axhline(y=power_limit, color='green', linestyle='-', linewidth=2.5, label='Limite (30kW)', zorder=5)

# Evidenzia violazioni PID
axs[1].fill_between(df['Time'], pid_power_motor, power_limit, 
                    where=(pid_power_motor > power_limit), 
                    color='red', alpha=0.25, label='Violazione PID', interpolate=True)

axs[1].set_xlabel('Tempo [s]', fontsize=13)
axs[1].set_ylabel('Potenza Motore [W]', fontsize=13)
axs[1].set_title('Gestione Vincolo Potenza', fontsize=15, fontweight='bold')
axs[1].legend(loc='upper right', fontsize=11, framealpha=0.9)
axs[1].grid(True, alpha=0.3)
axs[1].set_ylim([-2000, 35000])

plt.tight_layout()
plt.savefig('grafico_2_longitudinale.png', dpi=300, bbox_inches='tight')
print("Generato: grafico_2_longitudinale.png")

# ============================================================
# GRAFICO 3: FORZA LONGITUDINALE
# ============================================================
fig, ax = plt.subplots(figsize=(12, 5))

ax.plot(df['Time'], df['MPC_Fx'], 'b-', label='MPC Forza $F_x$', linewidth=2)
ax.plot(df['Time'], df['PID_Fx'], 'r-.', label='PID Forza $F_x$', linewidth=2)
ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax.axhline(y=1500, color='green', linestyle='--', linewidth=1.5, label='Forza Massima (1500N)', alpha=0.7)
ax.axhline(y=-1500, color='orange', linestyle='--', linewidth=1.5, label='Frenata Massima (-1500N)', alpha=0.7)

# Zone di frenata
ax.fill_between(df['Time'], df['MPC_Fx'], 0, where=(df['MPC_Fx'] < 0), color='blue', alpha=0.15, label='MPC Frenata')
ax.fill_between(df['Time'], df['PID_Fx'], 0, where=(df['PID_Fx'] < 0), color='red', alpha=0.15, label='PID Frenata')

ax.set_xlabel('Tempo [s]', fontsize=13)
ax.set_ylabel('Forza Longitudinale $F_x$ [N]', fontsize=13)
ax.set_title('Confronto Dinamica Longitudinale: Accelerazione e Frenata', fontsize=15, fontweight='bold')
ax.legend(loc='upper right', fontsize=11, framealpha=0.9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('grafico_3_forza_longitudinale.png', dpi=300, bbox_inches='tight')
print("Generato: grafico_3_forza_longitudinale.png")

plt.show()

# ============================================================
# STATISTICHE
# ============================================================
print("\n" + "="*60)
print("STATISTICHE COMPARAZIONE")
print("="*60)

# Errore laterale
mpc_lateral_error = np.sqrt(np.mean((df['MPC_Y'] - df['Ref_Y'])**2))
pid_lateral_error = np.sqrt(np.mean((df['PID_Y'] - df['Ref_Y'])**2))
print(f"Errore Laterale RMSE:")
print(f"  MPC: {mpc_lateral_error:.3f} m")
print(f"  PID: {pid_lateral_error:.3f} m")
print(f"  Miglioramento: {((pid_lateral_error - mpc_lateral_error) / pid_lateral_error * 100):.1f}%")

# Errore velocità
mpc_speed_error = np.sqrt(np.mean((df['MPC_Vx'] - df['Ref_Vx'])**2))
pid_speed_error = np.sqrt(np.mean((df['PID_Vx'] - df['Ref_Vx'])**2))
print(f"\nErrore Velocità RMSE:")
print(f"  MPC: {mpc_speed_error:.3f} m/s")
print(f"  PID: {pid_speed_error:.3f} m/s")
print(f"  Miglioramento: {((pid_speed_error - mpc_speed_error) / pid_speed_error * 100):.1f}%")

# Violazioni potenza
pid_violations = np.sum(pid_power_motor > power_limit)
pid_violation_time = pid_violations * 0.05
mpc_violations = np.sum(mpc_power_motor > power_limit)

print(f"\nViolazioni Vincolo Potenza (30kW):")
print(f"  MPC: {mpc_violations} steps")
print(f"  PID: {pid_violations} steps ({pid_violation_time:.2f}s)")

# Energia consumata
mpc_energy = np.trapz(mpc_power_motor, df['Time']) / 1000
pid_energy = np.trapz(pid_power_motor, df['Time']) / 1000
print(f"\nEnergia Consumata:")
print(f"  MPC: {mpc_energy:.1f} kJ")
print(f"  PID: {pid_energy:.1f} kJ")

print("="*60)
