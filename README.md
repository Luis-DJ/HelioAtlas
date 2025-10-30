🪐 HelioAtlas — An Interactive Journey of Comet 3I/ATLAS

HelioAtlas is a Python-based heliocentric simulator for the interstellar comet 3I/ATLAS.

It visualises the comet’s 2025–2026 passage through the inner Solar System using:

A live, top-down orbital view (Sun, Earth, Mars, and the comet)

Distance curves over time (Sun–comet, Earth–comet, Mars–comet)

Angular geometry (phase angle, elongation, ecliptic longitude/latitude)

An interactive timeline slider with contextual annotations and a day counter

<p align="center"> <img src="media/3I_ATLAS_Timeline.png" width="700" alt="3I/ATLAS Timeline"/> </p>
🌌 Why HelioAtlas?

HelioAtlas is both an educational and exploratory tool:

It lets you see where 3I/ATLAS is relative to Earth and Mars on any given day.

It shows how close the comet gets to the Sun (perihelion) and planets.

It computes and displays viewing geometry (phase and elongation), which are critical for understanding when/why a comet might be observable.

⚙️ Installation and Setup
1️⃣ Clone or download this repository

If you have Git installed:

git clone https://github.com/Luis-DJ/HelioAtlas.git
cd HelioAtlas


Or simply click the green Code button → Download ZIP, and extract it locally.

2️⃣ Install dependencies

Use the requirements.txt file to install all required packages:

pip install -r requirements.txt


This will install:

numpy, pandas, matplotlib

astroquery, astropy

ipywidgets

3️⃣ Run HelioAtlas

To launch the application:

python helioatlas.py


When launched, you’ll see:

A heliocentric orbit view (Sun, Earth, Mars, 3I/ATLAS)

A distance vs. time chart

A timeline slider with monthly markers and real-time readouts

🪐 Optional: Create a Windows executable

To share with others without Python:

pip install pyinstaller
pyinstaller --onefile helioatlas.py


This will create dist/helioatlas.exe, ready to run on Windows.

🛰 Data Sources

NASA JPL HORIZONS — accurate heliocentric positions (Sun, Earth, Mars, comet)

Minor Planet Center (MPC) — orbital elements for hyperbolic propagation

Wikimedia Commons — open-use imagery for educational visuals

🧮 Key Terms
Symbol	Description
λ (Lambda)	Heliocentric ecliptic longitude (deg)
β (Beta)	Heliocentric ecliptic latitude (deg)
Phase (S–C–E)	Sun–Comet–Earth angle (illumination geometry)
Elongation (S–E–C)	Sun–Earth–Comet angle (sky separation)
Perihelion	Closest approach to Sun — 2025-10-29, q ≈ 1.36 AU
Inclination	~175°, retrograde orbit
Eccentricity	6.142 (hyperbolic / interstellar)
🙌 Acknowledgments / Agradecimientos

English:
HelioAtlas was conceived and developed by Luis D. Jimenez, blending a lifelong passion for science, travel, and visual storytelling.
The project was created in collaboration with ChatGPT (OpenAI, GPT-5), which assisted with orbital modeling, data visualization, and documentation.

Special thanks to:

NASA JPL HORIZONS and Minor Planet Center (MPC) for public ephemeris data

Wikimedia Commons for open-license imagery used in educational materials

Español:
HelioAtlas fue concebido y desarrollado por Luis D. Jimenez, combinando una pasión de toda la vida por la ciencia, los viajes y la narración visual.
El proyecto fue creado en colaboración con ChatGPT (OpenAI, GPT-5), que brindó asistencia en el modelado orbital, la visualización de datos y la documentación.

Agradecimientos especiales a:

NASA JPL HORIZONS y Minor Planet Center (MPC) por los datos públicos de efemérides

Wikimedia Commons por las imágenes de libre uso empleadas con fines educativos

📜 License

This project is released under the MIT License.
See LICENSE
 for details.

You are free to use, modify, and share HelioAtlas for educational or commercial purposes, provided attribution to Luis D. Jimenez is maintained.

“Built by curiosity and collaboration — a celebration of human + AI creativity in astronomy.”