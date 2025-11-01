# 🪐 HelioAtlas — An Interactive Journey of Comet 3I/ATLAS

HelioAtlas is a Python-based heliocentric simulator for the interstellar comet **3I/ATLAS**.

It visualizes the comet’s 2025–2026 passage through the inner Solar System through:

- 🌀 A live, top-down orbital view (Sun, Earth, Mars, and the comet)
- 📈 Distance curves over time (Sun–Comet, Earth–Comet, Mars–Comet)
- 🌗 Angular geometry (phase angle, elongation, ecliptic longitude/latitude)
- 🕒 An interactive timeline slider with contextual annotations and a day counter

<p align="center">
  <img src="media/3I_ATLAS_Timeline.png" width="700" alt="3I/ATLAS Timeline"/>
</p>

---

## 🌌 Why HelioAtlas?

HelioAtlas is both an educational and exploratory tool:

- See where **3I/ATLAS** is relative to Earth and Mars on any given day.  
- Track how close the comet gets to the Sun (perihelion) and nearby planets.  
- Study viewing geometry — *phase* and *elongation* — to understand when/why the comet might be visible.

---

## ⚙️ Installation and Setup

### 1️⃣ Clone or download this repository

If you have Git installed:
```bash
git clone https://github.com/Luis-DJ/HelioAtlas.git
cd HelioAtlas
```

Or click the green **Code → Download ZIP** button and extract it locally.

---

### 2️⃣ Install dependencies

Use the provided `requirements.txt` file:

```bash
pip install -r requirements.txt
```

This installs:
- `numpy`, `pandas`, `matplotlib`
- `astroquery`, `astropy`
- `ipywidgets`

---

### 3️⃣ Run HelioAtlas

To launch the application:
```bash
python HelioAtlas.py
```

You’ll see:
- A heliocentric orbit view (Sun, Earth, Mars, Comet)
- A distance vs. time chart
- An interactive timeline slider with monthly markers and real-time readouts

---

## 💫 Download the Standalone App (Windows)

The latest compiled `.exe` version can be downloaded here:

➡️ **[HelioAtlas v0.1-preview Release](https://github.com/Luis-DJ/HelioAtlas/releases)**

**Notes**
- 🪟 No installation needed — just run `HelioAtlas_Comet.exe`.  
- ⚠️ Windows SmartScreen may warn since the EXE is unsigned. Click “More info → Run anyway.”  
- 🕐 First launch may take **10–60 seconds** on slower PCs (PyInstaller unpacks in the background).  

---

## 🧑‍💻 Building the EXE Yourself

If you’d like to compile your own version:

1. Install PyInstaller in your environment:
   ```bash
   pip install pyinstaller
   ```

2. Build the executable:
   ```bash
   pyinstaller HelioAtlas.py        --name HelioAtlas        --onefile        --windowed        --icon=HelioAtlas_icon.ico        --collect-data astroquery        --collect-data astropy
   ```

3. Your EXE will appear in:
   ```
   dist/HelioAtlas.exe
   ```

**Tip:**  
When packaging interactive Matplotlib apps, ensure you include this snippet at the top of your code:

```python
import sys, matplotlib
if getattr(sys, 'frozen', False):
    matplotlib.use("Qt5Agg")
else:
    try:
        matplotlib.use("TkAgg")
    except Exception:
        pass
```

This ensures GUI rendering works correctly in the compiled EXE.

---

## 🧩 Developer Tip — Windows Icon Cache

If your custom icon doesn’t appear after rebuilding the EXE:

- Rename the file (e.g., `HelioAtlas_v0.2.exe`), or  
- Clear the Windows icon cache:
  ```powershell
  ie4uinit.exe -ClearIconCache
  ```
  (or delete `%LocalAppData%\IconCache.db` and restart Explorer)

---

## 🛰 Data Sources

- **NASA JPL HORIZONS** — accurate heliocentric positions (Sun, Earth, Mars, Comet)  
- **Minor Planet Center (MPC)** — orbital elements for hyperbolic propagation  
- **Wikimedia Commons** — open-use imagery for educational visuals  

---

## 🧮 Key Terms

| Symbol | Description |
|:-------:|:------------|
| λ (Lambda) | Heliocentric ecliptic longitude (deg) |
| β (Beta) | Heliocentric ecliptic latitude (deg) |
| Phase (S–C–E) | Sun–Comet–Earth angle (illumination geometry) |
| Elongation (S–E–C) | Sun–Earth–Comet angle (sky separation) |
| Perihelion | Closest approach to Sun — 2025-10-29, q ≈ 1.36 AU |
| Inclination | ~175°, retrograde orbit |
| Eccentricity | 6.142 (hyperbolic / interstellar) |

---

## 🙌 Acknowledgments / Agradecimientos

**English:**  
HelioAtlas was conceived and developed by **Luis D. Jimenez**, blending a lifelong passion for science, travel, and visual storytelling.  
Created in collaboration with **ChatGPT (OpenAI, GPT-5)**, which assisted with orbital modeling, data visualization, and documentation.

Special thanks to:
- **NASA JPL HORIZONS** and **MPC** for open ephemeris data  
- **Wikimedia Commons** for open-license imagery

**Español:**  
HelioAtlas fue concebido y desarrollado por **Luis D. Jimenez**, combinando una pasión de toda la vida por la ciencia, los viajes y la narración visual.  
Creado en colaboración con **ChatGPT (OpenAI, GPT-5)**, que brindó asistencia en modelado orbital, visualización de datos y documentación.

Agradecimientos especiales a:
- **NASA JPL HORIZONS** y **Minor Planet Center (MPC)** por los datos públicos  
- **Wikimedia Commons** por imágenes de libre uso

---

## 📜 License

This project is released under the **MIT License**.  
See [`LICENSE`](LICENSE) for details.

You are free to use, modify, and share HelioAtlas for educational or commercial purposes, provided attribution to **Luis D. Jimenez** is maintained.

> “Built by curiosity and collaboration — a celebration of human + AI creativity in astronomy.”
