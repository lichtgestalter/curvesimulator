"""
Ich habe ein Programm geschrieben, dass die Orbitparameter der Koerper eines Sternensystems als Input erhaelt und daraus fuer
jeden Koerper die initiale Position und Geschwindigkeit (state vector) berechnet.

Mein Programm funktioniert gut für Sternensysteme mit einem dominierenden Körper, der viel mehr Masse hat als alle anderen im System.
Aber bei Systemen mit zwei oder mehr Sternen errechnet das Programm zu hohe Startgeschwindigkeiten, so dass die Sterne auseinander
fliegen statt umeinander zu kreisen.

Ich vermute, der Fehler liegt in meiner Berechnung des “state vector” oder des “gravitational parameter”.

Zumindest für Systeme mit genau zwei Körpern weiß ich, wie groß der Fehler ist und könnte ihn mit einem Korrekturfaktor bereinigen,
allerdings ohne zu ahnen, was ich da tue. Wenn das System aus zwei Sternen A und B besteht und A mehr Masse hat als B,
dann kann ich aus der zunächst berechneten Geschwindigkeit v mit dieser Formel die korrekte Geschwindigkeit ermitteln:

v / (1 + (Masse B / Masse A))

Was habe ich falsch gemacht? Wie mache ich es richtig?

Anbei als Auszug aus meinem Programm die Funktionen zur Berechnung des state vector und des gravitational parameter:
"""


def keplerian_elements_to_state_vector(self):
    """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
    Returns also true anomaly, eccentric anomaly, mean anomaly and the time of periapsis.
    [a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
    [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
    [c]: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors/19335#19335
    [d]: https://space.stackexchange.com/questions/55356/how-to-find-eccentric-anomaly-by-mean-anomaly
    [e]: https://github.com/alfonsogonzalez/AWP/blob/main/src/python_tools/numerical_tools.py
    Numbers in comments refer to numbered formulas in [a] and [b].
    Code based on [c]. Added calculation of eccentric anomaly based on the explanations
    in [d] using a stripped down version of [e]."""

    self.calc_orbit_angles()  # Ω, ω, ϖ
    self.calc_period_or_semi_major_axis()  # P, a
    self.calc_anomalies()  # L, ma, ea, nu, T
    P, a, e, i, Ω, ω, ϖ, L = self.P, self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
    ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas

    r = a * (1 - e * math.cos(ea))  # 4b: radius r
    h = math.sqrt(mu * a * (1 - e ** 2))  # 5b: specific angular momentum h
    x = r * (math.cos(Ω) * math.cos(ω + nu) - math.sin(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component x
    y = r * (math.sin(Ω) * math.cos(ω + nu) + math.cos(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component y
    z = r * (math.sin(i) * math.sin(ω + nu))  # 6b: position component z
    p = a * (1 - e ** 2)  # 7b: Semi-latus rectum. Used in velocity calculation.
    dx = (x * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.cos(Ω) * math.sin(ω + nu) + math.sin(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component x
    dy = (y * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.sin(Ω) * math.sin(ω + nu) - math.cos(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component y
    dz = (z * h * e / (r * p)) * math.sin(nu) + (h / r) * (math.cos(ω + nu) * math.sin(i))  # 7b: velocity component z
    return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors


def gravitational_parameter(bodies, g):
    """Calculate the gravitational parameter of masses orbiting a common barycenter
    https://en.wikipedia.org/wiki/Standard_gravitational_parameter"""
    mass = 0.0
    for body in bodies:
        mass += body.mass
    # print(f"Gravitational parameter {g * mass:.3f}")
    return g * mass

"""
Antwort perplexity:
Der Hauptgrund für die zu hohen Startgeschwindigkeiten in Systemen mit mehreren massereichen Körpern liegt in der Berechnung des gravitationalen Parameters. Ihre aktuelle Implementierung berücksichtigt nicht die relative Bewegung der Körper zueinander. Hier sind die Hauptprobleme und Lösungsansätze:
Gravitationaler Parameter
Die Funktion gravitational_parameter summiert einfach die Massen aller Körper auf. Dies ist nur für Systeme mit einem dominanten Zentralkörper korrekt. Für Systeme mit mehreren massereichen Körpern muss der effektive gravitationale Parameter für jeden Körper individuell berechnet werden.
Für ein Zwei-Körper-System sollte der gravitationale Parameter wie folgt berechnet werden:
python
def gravitational_parameter_two_body(m1, m2, g):
    return g * (m1 + m2)
Für Systeme mit mehr als zwei Körpern wird die Berechnung komplexer und erfordert möglicherweise numerische Methoden.
Berechnung des Zustandsvektors
Die Funktion keplerian_elements_to_state_vector geht von einem Zwei-Körper-Problem aus, bei dem ein Körper um einen fixen Zentralkörper kreist. Für Systeme mit mehreren massereichen Körpern muss dieser Ansatz angepasst werden.
Korrektur für Zwei-Körper-Systeme

Ihre Korrekturformel für Zwei-Körper-Systeme ist korrekt
 Diese Formel berücksichtigt die Bewegung beider Körper um den gemeinsamen Schwerpunkt.

Lösungsansätze
Individuelle Berechnung: Berechnen Sie für jeden Körper einen individuellen gravitationalen Parameter, der die Massen der anderen Körper berücksichtigt.
Schwerpunktbetrachtung: Berechnen Sie die Bewegung relativ zum Schwerpunkt des Systems statt zu einem fixen Zentralkörper.
N-Körper-Simulation: Für komplexere Systeme mit mehr als zwei Körpern sollten Sie eine vollständige N-Körper-Simulation in Betracht ziehen, die die gegenseitigen Gravitationswechselwirkungen aller Körper berücksichtigt.
Iterative Methoden: Verwenden Sie iterative Methoden zur Bestimmung stabiler Anfangsbedingungen für Mehr-Körper-Systeme.
Beispiel für eine verbesserte Implementierung
Hier ist ein Beispiel für eine verbesserte gravitational_parameter Funktion für ein Zwei-Körper-System:
python
def gravitational_parameter_improved(body, other_body, g):
    return g * (body.mass + other_body.mass)
Für die Berechnung des Zustandsvektors müssten Sie dann die relative Bewegung der Körper zueinander berücksichtigen. Dies erfordert eine umfassendere Überarbeitung Ihrer keplerian_elements_to_state_vector Funktion.
Beachten Sie, dass für Systeme mit mehr als zwei Körpern die Berechnung deutlich komplexer wird und möglicherweise numerische Methoden oder eine vollständige N-Körper-Simulation erforderlich sind12.
"""
