import math
import tkinter as tk
from tkinter import ttk, messagebox
import mysql.connector

def get_conn():
    return mysql.connector.connect(
        host="localhost",
        user="root",
        password="phdpython22!",
        database="brake_analysis"
    )

def save_simulation_to_db(material, params, summary):
    try:
        conn = get_conn()
        cursor = conn.cursor()

        sql = """
        INSERT INTO simulations
        (material, rho, k, alpha, surface_peak, back_peak, avg_final, peak_stress)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
        """

        vals = (
            material,
            params["rho"],
            params["k"],
            params["alpha"],
            summary["surface_peak"],
            summary["back_peak"],
            summary["avg_final"],
            summary["peak_stress"]
        )

        cursor.execute(sql, vals)
        conn.commit()
        cursor.close()
        conn.close()

    except Exception as e:
        print("Failed to save simulation:", e)


# PHYSICS FUNCTIONS

def braking_energy(M, v):
    return 0.5 * M * v**2


def disc_area_from_diameter(D_outer, inner_frac=0.5):
    Ro = D_outer / 2
    Ri = inner_frac * Ro
    return math.pi * (Ro**2 - Ri**2), Ro, Ri


def solve_through_thickness(
    e, rho, c, k, q_time_func,
    t_total, dt, nz, h_surface, h_back,
    T_init, T_amb
):
    dz = e / (nz - 1)
    alpha_th = k / (rho * c)
    stability_dt = 0.5 * dz**2 / alpha_th

    if dt > stability_dt:
        dt = stability_dt * 0.9

    Nt = int(math.ceil(t_total / dt)) + 1
    times = [i * dt for i in range(Nt)]
    times[-1] = t_total

    T = [T_init] * nz
    T_record = [T.copy()]
    r = alpha_th * dt / (dz**2)

    for it in range(1, Nt):
        t = times[it]
        T_new = T.copy()

        for i in range(1, nz - 1):
            T_new[i] = T[i] + r * (T[i + 1] - 2 * T[i] + T[i - 1])

        q_surf = q_time_func(max(0, t - dt))

        conduction_flux = k * (T[1] - T[0]) / dz
        T_new[0] = T[0] + dt / (rho * dz * c) * \
            (conduction_flux + q_surf - h_surface * (T[0] - T_amb))

        conduction_flux_back = k * (T[-2] - T[-1]) / dz
        T_new[-1] = T[-1] + dt / (rho * dz * c) * \
            (conduction_flux_back - h_back * (T[-1] - T_amb))

        T = T_new
        T_record.append(T.copy())

    return times, T_record


def calculate_thermal_stress_profile(T_profile, T_ref, alpha, E, nu=0.3):
    return [E * alpha * (T - T_ref) / (1 - nu) for T in T_profile]


class App:
    def __init__(self, root):
        self.root = root
        root.title("Brake Disc Thermal Analysis — Custom Materials")
        root.geometry("950x480")

        self.last_result = None
        self.build_ui()

    def build_ui(self):
        frm = ttk.LabelFrame(self.root, text="Input Parameters", padding=10)
        frm.pack(side="left", fill="y", padx=10, pady=10)

        self.entries = {}

        # Input parameters including material properties
        params = [

        ]

        for lbl, val in params:
            row = ttk.Frame(frm)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=lbl, width=22).pack(side="left")
            ent = ttk.Entry(row, width=12)
            ent.insert(0, val)
            ent.pack(side="right")
            self.entries[lbl] = ent

        ttk.Button(
            frm, text="Run Simulation",
            command=self.run_simulation
        ).pack(pady=10)

        # Results Frame
        res = ttk.LabelFrame(self.root, text="Results", padding=10)
        res.pack(side="right", fill="both", expand=True, padx=10, pady=10)

        self.results_text = tk.Text(res, wrap="word", width=60, height=34)
        self.results_text.pack(side="left", fill="both", expand=True)

        scroll = ttk.Scrollbar(res, orient="vertical",
                               command=self.results_text.yview)
        scroll.pack(side="right", fill="y")
        self.results_text.configure(yscrollcommand=scroll.set)

    def run_simulation(self):
        try:
            M = 5000
            v = 200
            D = 1
            t_brake = 60
            c = 460
            E = 2 * 1e9
            h = 100
            T_amb = 24

            rho = 30
            k = 20
            alpha = 10

            e = 0.05
            nz = 60
            t_cool = 60
            inner_frac = 0.5

            A_contact, Ro, Ri = disc_area_from_diameter(D, inner_frac)
            Ek = braking_energy(M, v)
            Q_disc = Ek / 2
            q0 = Q_disc / (A_contact * t_brake)

            # --- q_time defined BEFORE usage ---
            def q_time(t):
                return q0 if t <= t_brake else 0

            alpha_th = k / (rho * c)
            dz = e / (nz - 1)
            dt = min(0.05, 0.5 * dz**2 / alpha_th * 0.9)

            times, T_record = solve_through_thickness(
                e, rho, c, k,
                q_time,
                t_brake + t_cool,
                dt,
                nz,
                h, h,
                T_amb, T_amb
            )

            final_T = T_record[-1]
            surface_peak = max(t[0] for t in T_record)
            back_peak = max(t[-1] for t in T_record)
            avg_final = sum(final_T) / nz
            stress = calculate_thermal_stress_profile(final_T, T_amb, alpha, E)
            peak_stress = max(abs(s) for s in stress) / 1e6

            params_dict = {"rho": rho, "k": k, "alpha": alpha}
            summary = {
                "surface_peak": surface_peak,
                "back_peak": back_peak,
                "avg_final": avg_final,
                "peak_stress": peak_stress
            }

            # Save into DB as "Custom"
            save_simulation_to_db("Custom", params_dict, summary)

            # Display result
            txt = f"""
=== SIMULATION RESULTS ===

ρ: {rho}
k: {k}
α: {alpha}

Peak surface T: {surface_peak:.2f} °C
Peak back-face T: {back_peak:.2f} °C
Final average T: {avg_final:.2f} °C

Peak thermal stress: {peak_stress:.2f} MPa

(Simulation saved to MySQL)
"""
            self.results_text.delete("1.0", "end")
            self.results_text.insert("1.0", txt)

        except Exception as e:
            messagebox.showerror("Error", f"Input error:\n{e}")

if __name__ == "__main__":
    root = tk.Tk()
    App(root)
    root.mainloop()