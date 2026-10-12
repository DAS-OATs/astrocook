# Astrocook – Roadmap di sviluppo

Documento di lavoro, nato da un'analisi del codice v2.0.1 (ottobre 2026).
Raccoglie la valutazione dello stato attuale, il piano di consolidamento
della v2 e la visione per una futura v3. Va aggiornato alla fine di ogni fase.

---

## 1. Stato attuale (v2.0.1)

### Cos'è
Pacchetto Python per l'analisi di spettri di quasar, con GUI PySide6 e API
scriptabile. Copre I/O (FITS, ESO/MIDAS, archivi `.acs` V1 e `.acs2` V2),
operazioni sul flusso, stima del continuo, identificazione delle righe e
fit di Voigt (minimi quadrati + motore bayesiano sperimentale emcee/dynesty).

Dimensioni: ~53k righe Python, di cui ~26k di codice V2 attivo e ~22.5k in
`astrocook/legacy/`.

### Struttura
- `astrocook/core/` – strutture dati immutabili (`structures.py`),
  `SpectrumV2` (`spectrum.py`), funzioni pure (`spectrum_operations.py`),
  sessione e lista dei sistemi.
- `astrocook/fitting/` – `voigt_fitter.py`, `voigt_model.py` (vincoli),
  `bayesian_fitter.py`.
- `astrocook/recipes/` – operazioni esposte all'utente (continuum,
  absorbers, edit, file, flux), con parametri strutturati.
- `astrocook/gui/` – interfaccia PySide6.
- `astrocook/io/` – loader e migrazione V1 → V2.
- `astrocook/legacy/` – codice V1 (wxPython), ancora in parte importato.

### Punti di forza
1. Strutture dati immutabili (dataclass frozen): undo affidabile, cronologia
   coerente, niente effetti collaterali.
2. Unità fisiche esplicite (astropy) su ogni colonna.
3. Separazione a livelli: core/fitting, recipe, GUI.
4. Riproducibilità: log di sessione con scripter; migrazione automatica V1.
5. Fit di Voigt maturo, con motore bayesiano e incertezze asimmetriche.
6. Docstring numpydoc, documentazione Sphinx con tutorial, CHANGELOG curato,
   build macOS e Docker.

### Punti di debolezza
1. **Legacy ancora accoppiato:** 17 moduli V2 (anche `core/session.py` e
   `core/spectrum.py`) importano da `astrocook.legacy` (`Defaults`, `Format`,
   `GUILog`, `SpectrumV1`, `SystListV1`, `trans_parse`, `xem_d`,
   `qso_composite`, …).
2. **Test scarsi:** ~5 funzioni `test_*`; il resto di `tests/` sono script di
   verifica. Nessuna CI. Nessun test di regressione numerica.
3. **File molto grandi:** `gui/main_window.py` (4.6k righe),
   `gui/pyside_plot.py` (2.9k), `core/spectrum.py` (2.8k),
   `gui/system_inspector.py` (2.1k). `SpectrumV2` è una "god class".
4. **Gestione degli errori larga:** ~200 `except:` / `except Exception`
   fuori dal legacy. API incoerenti (es. `load_session_from_file` restituisce
   `0` in caso di errore).
5. **Packaging:** manca il file LICENSE (campo commentato in
   `pyproject.toml`); placeholder del template ancora presenti; entry point
   `astrocook-gui` commentato; `numpy>=1.17.3` non realistico con
   Python ≥3.12; emcee/dynesty/corner obbligatori invece che opzionali.
6. **Disordine nel repository:** `docs_backup/`, file `*_old.py`, `*.bck`,
   `telluric_old.fits`, `atom_par_backup.dat` (alcuni inclusi nel pacchetto).

---

## 2. Valutazione scientifica

### Fit di Voigt – solido
- Profilo esatto via funzione di Faddeeva (`wofz`); b termico e turbolento
  in quadratura.
- Convoluzione gaussiana con R costante, variabile per pixel o per componente.
- `least_squares` (TRF) con limiti fisici e `x_scale='jac'`; vincoli tra
  componenti; aggiunta iterativa di componenti con criterio AIC.
- Motore bayesiano (emcee, dynesty con evidenza).

**Da verificare / mancante:**
1. **Nessun sovracampionamento:** il modello è calcolato al centro dei pixel.
   Per righe strette rispetto al pixel può distorcere b e N.
   *Primo controllo da fare.*
2. **Incertezze non scalate per il χ² ridotto:** il commento in
   `fitting/voigt_fitter.py` (~riga 936) dice di farlo, il codice no.
   Da chiarire se è voluto. *Primo controllo da fare.*
3. Solo profilo strumentale gaussiano (manca un LSF tabulato, utile per
   ESPRESSO/UVES).
4. Mancano offset del livello zero e continuo locale libero nella finestra
   di fit.

### Continuo – punto debole
Algoritmo attuale (`find_absorbed_regions` + `fit_continuum`): media mobile
con sigma-clipping iterativo (scale diverse dentro/fuori dalla foresta Lyα),
interpolazione lineare sulle regioni assorbite, fattore "fudge" globale,
smoothing gaussiano.

Limiti: nella foresta ad alto z il livello parte troppo basso e un fattore
globale non corregge un errore che cresce con z; euristica `/8` ereditata da
V1 per la finestra; template QSO solo composito, PCA non implementata;
nessuna incertezza sul continuo né propagazione a N, b.

Direzioni (sforzo crescente):
1. Continuo locale (polinomio di grado basso) fittato insieme alle righe.
2. Spline con nodi adattivi e incertezza (es. processo gaussiano).
3. Predizione del continuo della foresta dal lato rosso (PCA/ML su mock).
4. Correzione basata sull'evoluzione dell'assorbimento medio con z.

### Identificazione dei sistemi – sensata ma limitata
Algoritmo attuale: prodotto delle significatività su griglia di z comune
(`compute_identification_signal`, da V1); per i doppietti, sovrapposizione
cinematica + punteggio sulla relazione F₁ = F₂^γ (`rate_doublet_candidate`).

Limiti: solo le prime due righe di ogni multipletto, una serie alla volta;
nessuna risoluzione globale dei conflitti; F₁ = F₂^γ vale per il profilo
intrinseco, non dopo la convoluzione; esito sì/no, non probabilistico;
contaminanti Lyα non modellati.

Direzioni:
1. Identificazione probabilistica come confronto tra ipotesi.
2. Ciclo identifica → fitta → residui → identifica.
3. Rivelatori addestrati su mock per specie specifiche (CIV, MgII, DLA).

### Validazione
Prima di nuovi algoritmi: validazione quantitativa con spettri simulati
(parametri noti, vari S/N e risoluzioni) e confronto con VPFIT o altri
strumenti. Possibile base per un articolo su Astrocook 2.

---

## 3. Consolidamento della v2

**Principio:** prima si costruisce la rete di sicurezza, poi si tocca il
codice. Ogni passo è una PR piccola, con la suite verde prima e dopo.

### Fase 0 – Ambiente
Installazione delle dipendenze; avvio della GUI con
`QT_QPA_PLATFORM=offscreen`.

### Fase 1 – Rete di sicurezza (pilota)
Nessuna modifica al codice di produzione.
1. **Snapshot numerici** (characterization test) delle operazioni di
   `SpectrumV2` e delle funzioni pure: output del codice attuale salvati e
   confrontati con tolleranze esplicite.
2. **Replay dei log di sessione** come test end-to-end (servono sessioni
   reali fornite da Guido).
3. **Spettri simulati** con parametri noti (estendere
   `tests/simulate_and_fit_cii.py`). Progettare il generatore pensando anche
   a completezza/purezza per la v3.
4. **Non determinismo:** seed fissi, tolleranze statistiche, esclusione di
   UUID e timestamp.
5. **Compatibilità file:** round-trip `.acs`/`.acs2` (servono archivi V1
   reali).
6. **Misura della rete:** coverage; mutation testing (`mutmut`) sui moduli
   critici.
7. **GUI:** smoke test con `pytest-qt` in modalità offscreen.
8. **CI:** GitHub Actions (Python 3.12, 3.13).
9. **Primi controlli scientifici:** sovracampionamento del modello di Voigt;
   scalatura delle incertezze per il χ² ridotto.

A fine fase: misurare consumo e tempo di revisione per ristimare il resto.

### Fase 2 – Pulizia e packaging
File obsoleti e backup; `pyproject.toml` (placeholder, vincoli, entry point,
extra opzionali `[bayes]`); **LICENSE (decisione di Guido)**.

### Fase 3 – Distacco dal legacy
Una dipendenza alla volta (`Defaults`, `trans_parse`, `xem_d`, … diventano
moduli V2). Il legacy resta solo dietro il loader di migrazione V1.

### Fase 4 – Suddivisione dei file grandi
`SpectrumV2` in moduli tematici (continuo, righe, flusso) mantenendo la
stessa API pubblica tramite delega; stesso approccio per la GUI.

### Fase 5 – Gestione degli errori
Le `except` generiche, caso per caso: eccezione specifica, log esplicito,
oppure punto da discutere.

### Stime indicative
10–20 sessioni di lavoro, 2–5 settimane di calendario; il collo di bottiglia
è la revisione umana e la disponibilità di dati di test.

### Decisioni riservate a Guido
- Licenza.
- Comportamenti dubbi trovati dai characterization test (bug o voluti?).
- Validazione scientifica dei risultati.
- Release (PyPI, firma app macOS).

---

## 4. Direzione: Astrocook per gli agenti

Non integrare un LLM dentro Astrocook, ma rendere Astrocook facile da usare
per qualunque LLM.

1. **Server MCP:** ogni recipe come strumento con nome, descrizione e
   parametri tipizzati (le recipe hanno già parametri strutturati).
2. **API Python** pulita, documentata, con output strutturati (JSON:
   parametri con errori, χ², flag di qualità).
3. **Log affidabile in linguaggio naturale**, pronto da copiare in un testo
   (es. sezione metodi di un articolo). Registrare anche le *decisioni* del
   fit (componenti aggiunte/scartate e criterio), non solo i risultati.
4. **GUI come banco di verifica:** l'agente fa il lavoro ripetitivo, l'umano
   ispeziona, corregge, approva.

Pilota proposto: server MCP con 3–4 recipe, dopo la fase 1.

---

## 5. Visione v3 (orizzonte, dopo il consolidamento della v2)

**Obiettivo:** dato uno spettro, produrre in meno di un'ora **una collezione
di modelli completi**, con incertezze e una mappa di affidabilità, che
riportino anche le possibili diverse identificazioni della stessa feature,
ciascuna con un parametro che ne quantifichi la credibilità.

### Ciclo iterativo (deterministico e riproducibile)
1. Continuo iniziale con incertezza.
2. Rilevamento delle regioni assorbite e scomposizione in blocchi
   indipendenti (regioni collegate da multipletti = un blocco).
3. Identificazione come test di ipotesi: metalli, poi serie di Lyman,
   Lyα come ipotesi di riserva.
4. Fit di Voigt in parallelo per blocco, con continuo locale libero e
   aggiunta di componenti guidata da residui e criterio informativo.
5. Raffinamento del continuo usando il modello.
6. Ripetizione fino a convergenza.
7. Output: collezione di modelli, componenti con incertezze, mappa di
   affidabilità (χ² locale, alternative scartate e margine, allarmi).

### Credibilità dei modelli
- Pesi di Akaike/Schwarz (AIC/BIC) come stima rapida; evidenza bayesiana
  (dynesty) per i casi importanti.
- Alternative definite **per blocco**, non su tutto lo spettro, per evitare
  l'esplosione combinatoria.

### Difficoltà note
χ² ≈ 1 è necessario ma non sufficiente: degenerazione continuo/assorbimento
nella foresta; identificazioni sbagliate con χ² perfetto; sovra-adattamento;
spettro d'errore non affidabile (pixel correlati); casi patologici (DLA,
BAL, tellurici, artefatti).

### Validazione
Mock con verità nota: completezza, purezza, bias su N/b/z in funzione di
S/N e risoluzione, tasso di identificazioni errate.

### Percorso
Pilota sui metalli fuori dalla foresta, poi estensione alla foresta Lyα.
L'agente LLM resta ai margini (strategia nei casi bloccati, guida alla
revisione); ogni sua decisione finisce nel log.

### Implicazioni per il lavoro sulla v2
- Il generatore di mock della fase 1 è il primo nucleo del banco di prova v3.
- Nei refactoring non introdurre l'ipotesi "un solo modello per sessione".
- Il log deve registrare le decisioni del fit.

---

## 6. Stato di avanzamento

| Fase | Stato |
|---|---|
| Analisi iniziale e roadmap | completata (ottobre 2026) |
| 0 – Ambiente | da iniziare |
| 1 – Rete di sicurezza | da iniziare |
| 2 – Pulizia e packaging | da iniziare |
| 3 – Distacco dal legacy | da iniziare |
| 4 – Suddivisione file grandi | da iniziare |
| 5 – Gestione degli errori | da iniziare |
