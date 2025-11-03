# Substructure Search Service

Simple FastAPI application for storing molecules, running RDKit-powered
substructure searches, and exposing the results over a REST API.

## Features
- CRUD API for molecules identified by SMILES.
- RDKit substructure and similarity search utilities.
- Async search pipeline using Celery workers.
- Redis-backed caching to speed up repeated queries.
- Docker Compose stack with optional Nginx load balancer.
- Pytest suite and flake8 linting via GitHub Actions.

## Quick Start
1. Install Docker and Docker Compose.
2. Copy `.env.example` to `.env` and adjust connection settings if needed.
3. Build and start the stack:
   ```bash
   docker-compose up --build -d
   ```
4. Open `http://localhost:8000/docs` for interactive API docs.
5. Hit `http://localhost/` (port 80) to see load-balancer routing between web nodes.

## Local Development
```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
uvicorn src.main:app --reload
```

## Tests & Linting
```bash
pytest
flake8 src tests
```

## Useful Endpoints
- `POST /molecules` – add a molecule.
- `GET /molecules/{id}` – retrieve stored molecule.
- `GET /molecules?limit=100` – stream molecules with an iterator.
- `POST /search/substructure` – synchronous RDKit search.
- `POST /tasks/substructure` + `GET /tasks/{task_id}` – async search via Celery.

## Project Structure
- `src/` – FastAPI app, RDKit logic, Celery workers.
- `tests/` – unit and integration tests.
- `docker-compose.yml` – FastAPI, Redis, Celery, and optional Nginx services.
- `nginx/` – load-balancer configuration.
