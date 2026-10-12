# CLAUDE.md

Contesto per le sessioni di Claude Code su Astrocook. Il piano completo
(analisi, fasi, valutazione scientifica, visione v3) è in `ROADMAP.md`:
leggilo prima di iniziare un lavoro non banale.

## Progetto
Astrocook: pacchetto Python per l'analisi di spettri di quasar (continuo,
identificazione dei sistemi di assorbimento, fit di Voigt). GUI PySide6 +
API scriptabile. Autore e maintainer: Guido Cupani (INAF–OATs).
Versione attuale: 2.0.1. Python ≥3.12.

## Architettura (v2)
- `astrocook/core/` – dataclass immutabili (`structures.py`), `SpectrumV2`,
  funzioni pure in `spectrum_operations.py`, sessione, lista dei sistemi.
  I metodi restituiscono nuovi oggetti: non mutare lo stato in place.
- `astrocook/fitting/` – fit di Voigt (`voigt_fitter.py`), vincoli
  (`voigt_model.py`), motore bayesiano (`bayesian_fitter.py`).
- `astrocook/recipes/` – operazioni esposte all'utente, con parametri
  strutturati.
- `astrocook/gui/` – interfaccia PySide6.
- `astrocook/io/` – loader e migrazione V1 → V2.
- `astrocook/legacy/` – codice V1, ancora importato da vari moduli V2.
  Va rimosso progressivamente (fase 3), non usato per codice nuovo.

## Come si lavora
- Avvio GUI: `python run_astrocook.py [file ...]`
  (headless: `QT_QPA_PLATFORM=offscreen`).
- Test: `pytest` (config in `pyproject.toml`). La copertura è ancora
  minima: prima di modificare codice di produzione, aggiungi test di
  caratterizzazione per il comportamento coinvolto.
- Lavora per PR piccole; la suite deve essere verde prima e dopo.
- Non cambiare comportamenti numerici senza l'accordo di Guido: se un
  comportamento sembra un bug, segnalalo invece di correggerlo in silenzio.
- La documentazione utente (Sphinx + MyST) è in `docs/` ed è in inglese.

## Direzione
- v2: consolidamento (test e CI, pulizia, distacco dal legacy, refactoring,
  gestione degli errori).
- Astrocook come strumento per agenti: server MCP sulle recipe, API con
  output strutturati, log affidabile in linguaggio naturale.
- v3 (orizzonte): collezione di modelli completi con identificazioni
  alternative e credibilità, generata automaticamente in meno di un'ora.

## Stato attuale
Analisi e roadmap completate (ottobre 2026). Prossimo passo: fase 1 (rete di
sicurezza), come pilota. Aggiorna questa sezione e la tabella in
`ROADMAP.md` alla fine di ogni fase.
