# Code- und Mathe-Review (BGSC)

Dieses Dokument fasst die zwei Review-Runden zusammen (1. allgemeine
Code-Sauberkeit, 2. mathematische Korrektheit der Bayes'schen
Genklassifikation und der qPCR-ΔΔCt-Methode), und listet für jeden Punkt auf,
ob und wie er in Version 0.0.1.4 behoben wurde. Ziel war, echte Fehler zu
beheben und die Robustheit zu verbessern, **ohne** die publizierten
Zahlen/Ergebnisse der zugehörigen BMC-Bioinformatics-Publikation zu verändern.

## Zusammenfassung des Ergebnisses

Der mathematische Kern (Log-Likelihood → BIC → Posterior-Klassifikation,
log2FC ± Standardfehler für Microarray und qPCR) war bereits korrekt und
intern konsistent implementiert. Die gefundenen Probleme lagen entweder in
totem Code, in Diagnose-Nebenspalten oder in nicht genutzten Codepfaden – sie
haben die publizierten Ergebnisse nicht verändert, waren aber echte Bugs bzw.
Robustheitslücken, die jetzt behoben sind.

## Behobene Probleme

| # | Problem | Quelle | Datei | Status |
|---|---------|--------|-------|--------|
| 1 | `make.plot.data.exp.Ill.qPCR()`: `qgenesIDs` fehlte als Parameter, Funktion griff auf ein globales Objekt zu, das bei einer Paketfunktion **nicht** im Scope liegt → Aufruf in `ReproducibleScript.Rmd`/`scripts/MAIN.R` hätte mit "object 'qgenesIDs' not found" abgebrochen | Antwort 1 (Clean Code) | `R/BGSC.R` | ✅ behoben: Parameter ergänzt, beide Aufrufstellen + `man/*.Rd` angepasst |
| 2 | `get.Posterior()`: unnötiger exp→log→exp-Rundtrip ohne Log-Sum-Exp-Stabilisierung, Overflow-Risiko bei extremen BIC-Werten | Antwort 2 (Mathe) | `R/BGSC.R` | ✅ behoben: direkte Log-Space-Berechnung mit Log-Sum-Exp-Trick, Ergebnis für normale Wertebereiche identisch |
| 3 | `InformationCriterion(..., IC='AIC')`: `k` (Anzahl Beobachtungen) wurde fälschlich als AIC-Strafterm-Multiplikator verwendet statt der festen 2 | Antwort 2 (Mathe) | `R/BGSC.R` | ✅ behoben: `AIC = -2*logL + 2*npar`. **Kein Effekt auf publizierte Ergebnisse**, da im Pipeline-Code nur `IC='BIC'` verwendet wird (unverändert) |
| 4 | `InformationCriterion()`: Fehlerfälle nur `print()`-t statt echten Fehler zu werfen, Rückgabe `NULL` | Antwort 1 (Clean Code) | `R/BGSC.R` | ✅ behoben: `stop()` mit klarer Fehlermeldung |
| 5 | `comparativeMethod_qPCR.analysis()`: "Henri"-Zweig (`ddCTPower2`/`ddCTRelExp`) berechnete `2^(+ddCt)` statt der Standardkonvention `2^(-ddCt)` – Kehrwert der korrekten relativen Expression | Antwort 2 (Mathe) | `R/BGSC.R` | ✅ behoben: Vorzeichen korrigiert, stimmt jetzt mit dem unabhängig berechneten "Johanna"-Zweig (`Potenz`/`relative_values`) überein. **Kein Effekt auf publizierte Ergebnisse**, da diese Spalten nirgends weiterverwendet wurden (toter Code) |
| 6 | `comparativeMethod_qPCR.RNAi.log2FC()`: doppelt vergebener Spaltenname in `dCT.C1C0`/`relVal.C1C0` (C0-Standardfehler-Spalte hieß wie die C1-Spalte) | Antwort 2 (Mathe) | `R/BGSC.R` | ✅ behoben: Spalten korrekt benannt. **Kein Effekt auf publizierte Ergebnisse**, da `dCT.logFC.C1C0.stderr` bereits vorher aus den korrekten Variablen berechnet wurde |
| 7 | `makeExpData()`/`.makeQPCR()`: hartkodierter absoluter Pfad (`/Users/weinhol/GitHub/BGSC/`) und `setwd()`-Aufrufe | Antwort 1 (Clean Code) | `R/BGSC.R` | ✅ behoben: `pkgDir`-Parameter (Default `getwd()`), kein `setwd()` mehr. Betrifft nur interne, nicht exportierte Einmal-Hilfsfunktionen zur Datengenerierung, nicht die Analyse-Pipeline |
| 8 | Tote/auskommentierte Debug-Ausgaben (`# print(...)`), nie ausgeführter `FYI`-Block, ungenutzte Berechnung (`normexp.fit.detection.p`) | Antwort 1 (Clean Code) | `R/BGSC.R` | ✅ entfernt. Die im `FYI`-Block geprüfte Eigenschaft (Satterthwaite-Approximation == gepoolter Standardfehler) ist jetzt ein echter Test (`tests/testthat/test-log2FC.R`) |
| 9 | Keine Tests trotz `testthat` in `Suggests` | Antwort 1 (Clean Code) | – | ✅ `tests/testthat/` mit 6 Testdateien ergänzt (siehe unten) |
| 10 | **Neu gefunden** (nicht Teil der ursprünglichen zwei Antworten): `union5()` rief intern `Reduce(intersect, ...)` statt `Reduce(union, ...)` auf (Copy-Paste-Fehler von `intersect5()`) und lieferte damit die Schnittmenge statt der Vereinigung | – | `R/VennFunction.R` | ✅ behoben, da beim Umsetzen der Fixes entdeckt. Kein bekannter Aufrufer im Paket (exportierte, aber ungenutzte Funktion) |

## Bewusst nicht behoben (Begründung)

Diese Punkte aus Antwort 1 wurden **nicht** angefasst, um das Risiko zu
vermeiden, die Reproduzierbarkeit der publizierten Ergebnisse zu gefährden,
ohne die Möglichkeit zu haben, die Pipeline tatsächlich mit R auszuführen und
gegen die Originalzahlen zu verifizieren (in dieser Umgebung ist kein R
installiert):

- **Code-Duplikation** zwischen `comparativeMethod_qPCR.RNAi.log2FC()` und
  `comparativeMethod_qPCR.EGF.log2FC()` (~120 Zeilen sich stark ähnelnder
  Logik) – ein Refactoring in gemeinsame Hilfsfunktionen wäre sinnvoll, ist
  aber ein größerer strukturneller Eingriff.
- **Inkonsistente Namenskonventionen** (`camelCase` vs. `dot.case` vs.
  `snake_case`) über die gesamte exportierte API – Umbenennungen würden die
  öffentliche API brechen (`NAMESPACE`, bestehende Aufrufer/Skripte externer
  Nutzer der publizierten Version).
- **Magic Numbers** (`1:6`, `rep(0,6)` etc. für die feste Anzahl von 6 Proben)
  – zentrale Konstante wäre sauberer, ändert aber viele Stellen gleichzeitig.
- **`T`/`F` statt `TRUE`/`FALSE`** – rein stilistisch, kein Bug.
- **Reshape-Annahme** in `comparativeMethod_qPCR.EGF.log2FC()`
  (`matrix(x$dCT, nREP, nrow(x)/nREP)`): verlässt sich auf eine implizite
  Zeilenreihenfolge der qPCR-Rohdaten statt gegen die vorhandenen
  `Treatment`/`CellLine`-Spalten zu validieren. Ohne Zugriff auf die realen
  Rohdaten (nur `.rda`, kein R in dieser Umgebung verfügbar) konnte ich nicht
  verifizieren, ob/wie das zu beheben wäre, ohne das Verhalten für die
  tatsächlichen Daten zu verändern. Als Risiko dokumentiert, nicht als
  bestätigter Fehler.

## Neue Tests

- `tests/testthat/test-information-criterion.R` – BIC-Formel, feste
  AIC-Strafe (unabhängig von `k`), Fehlerfälle.
- `tests/testthat/test-posterior.R` – Bayes-Formel von `get.Posterior()`,
  Stabilität bei extremen BIC-Werten, Grenzfall gleicher BIC-Werte.
- `tests/testthat/test-loglikelihood.R` – `logLikelihoodOfnormData()` gegen
  manuell berechnete Gauß-Log-Likelihood, Freiheitsgrade der gepoolten
  Varianz.
- `tests/testthat/test-qpcr.R` – Vorzeichen-Regressionstest für den
  "Henri"-Fix, Kalibrator-Wert = 100.
- `tests/testthat/test-log2FC.R` – Satterthwaite-Äquivalenz,
  Spaltennamen-Regressionstest.
- `tests/testthat/test-venn-sets.R` – `intersect3/4/5`, `union3/4/5`
  (inkl. Regressionstest für den `union5`-Fix).

**Hinweis:** In dieser Umgebung ist kein R/`Rscript` installiert, die Tests
konnten daher nicht tatsächlich ausgeführt werden. Alle Formeln wurden von
Hand nachgerechnet (siehe vorherige Analyse-Antworten) und die R-Syntax
sorgfältig manuell geprüft; ein Lauf mit `devtools::test()` bzw.
`R CMD check` vor einem Merge wird dringend empfohlen.
