from fastapi import FastAPI, HTTPException, status, Depends, Query
from contextlib import asynccontextmanager
from pydantic import BaseModel
from os import getenv
from sqlalchemy.orm import Session
from sqlalchemy import inspect
from typing import Iterator, Optional
from rdkit import Chem
from celery.result import AsyncResult
from src.database import engine, get_db
from src.models import Base, Molecule
from src.logging_config import logger
from src.cache import get_cached_result
from src.celery_app import celery_app
from src.tasks import run_substructure_search



class MoleculeCreate(BaseModel):
    id: str
    smiles: str


class MoleculeUpdate(BaseModel):
    smiles: str


class MoleculeOut(BaseModel):
    id: str
    smiles: str


class SubstructureRequest(BaseModel):
    substructure: str


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("Starting application...")

    inspector = inspect(engine)
    existing_tables = inspector.get_table_names()

    if "molecules" not in existing_tables:
        try:
            Base.metadata.create_all(bind=engine)
            logger.info("Database tables created successfully")
        except Exception as e:
            inspector = inspect(engine)
            if "molecules" not in inspector.get_table_names():
                logger.error(f"Failed to create database tables: {e}")
                raise
            logger.info("Database tables already exist (created by another container)")
    else:
        logger.info("Database tables already exist")

    yield

    logger.info("Shutting down application...")

app = FastAPI(title="Molecule", lifespan=lifespan)


@app.get("/")
def get_server():
    """Endpoint to check load balancing - returns the server ID."""
    server_id = getenv("SERVER_ID", "1")
    logger.debug(f"Load balancer check - Server ID: {server_id}")
    return {"server_id": server_id}


def _validate_smiles(smiles: str) -> None:
    if Chem.MolFromSmiles(smiles) is None:
        logger.warning(f"Invalid SMILES provided: {smiles!r}")
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_CONTENT, detail=f"Invalid SMILES: {smiles!r}")




# Add molecule
@app.post("/molecules", response_model=MoleculeOut, status_code=status.HTTP_201_CREATED)
def add_molecule(m: MoleculeCreate, db: Session = Depends(get_db)):
    logger.info(f"POST /molecules - Adding molecule with ID: {m.id}")
    # Check if molecule with this ID already exists
    existing = db.query(Molecule).filter(Molecule.id == m.id).first()
    if existing:
        logger.warning(f"POST /molecules - Conflict: ID {m.id} already exists")
        raise HTTPException(status_code=status.HTTP_409_CONFLICT, detail="Identifier already exists")
    _validate_smiles(m.smiles)
    db_molecule = Molecule(id=m.id, smiles=m.smiles)
    db.add(db_molecule)
    db.commit()
    db.refresh(db_molecule)
    logger.info(f"POST /molecules - Successfully added molecule with ID: {m.id}")
    return MoleculeOut(id=db_molecule.id, smiles=db_molecule.smiles)


# Get molecule by id
@app.get("/molecules/{mol_id}", response_model=MoleculeOut)
def get_molecule(mol_id: str, db: Session = Depends(get_db)):
    logger.info(f"GET /molecules/{mol_id} - Retrieving molecule")
    db_molecule = db.query(Molecule).filter(Molecule.id == mol_id).first()
    if not db_molecule:
        logger.warning(f"GET /molecules/{mol_id} - Not found")
        raise HTTPException(status_code=404, detail="Not found")
    logger.debug(f"GET /molecules/{mol_id} - Successfully retrieved")
    return MoleculeOut(id=db_molecule.id, smiles=db_molecule.smiles)


# Updating a molecule by id
@app.put("/molecules/{mol_id}", response_model=MoleculeOut)
def update_molecule(mol_id: str, upd: MoleculeUpdate, db: Session = Depends(get_db)):
    logger.info(f"PUT /molecules/{mol_id} - Updating molecule")
    db_molecule = db.query(Molecule).filter(Molecule.id == mol_id).first()
    if not db_molecule:
        logger.warning(f"PUT /molecules/{mol_id} - Not found")
        raise HTTPException(status_code=404, detail="Not found")
    _validate_smiles(upd.smiles)
    db_molecule.smiles = upd.smiles
    db.commit()
    db.refresh(db_molecule)
    logger.info(f"PUT /molecules/{mol_id} - Successfully updated")
    return MoleculeOut(id=db_molecule.id, smiles=db_molecule.smiles)


# Delete a molecule by identifier
@app.delete("/molecules/{mol_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_molecule(mol_id: str, db: Session = Depends(get_db)):
    logger.info(f"DELETE /molecules/{mol_id} - Deleting molecule")
    db_molecule = db.query(Molecule).filter(Molecule.id == mol_id).first()
    if not db_molecule:
        logger.warning(f"DELETE /molecules/{mol_id} - Not found")
        raise HTTPException(status_code=404, detail="Not found")
    db.delete(db_molecule)
    db.commit()
    logger.info(f"DELETE /molecules/{mol_id} - Successfully deleted")
    return None


# List all molecules
@app.get("/molecules", response_model=list[MoleculeOut])
def list_molecules(
    limit: Optional[int] = Query(None, ge=1, description="Maximum number of molecules to return"),
    db: Session = Depends(get_db)
):
    logger.info(f"GET /molecules - Listing molecules (limit: {limit})")
    
    query = db.query(Molecule)
    if limit is not None:
        query = query.limit(limit)
    
    def molecule_iterator() -> Iterator[MoleculeOut]:
        for molecule in query:
            yield MoleculeOut(id=molecule.id, smiles=molecule.smiles)
    
    molecules = list(molecule_iterator())
    logger.debug(f"GET /molecules - Returning {len(molecules)} molecules")
    return molecules


@app.post("/search/substructure")
def start_substructure_search(payload: SubstructureRequest):
    substructure = payload.substructure
    logger.info("POST /search/substructure - Starting task for substructure: %s", substructure)

    if Chem.MolFromSmiles(substructure) is None:
        logger.warning("POST /search/substructure - Invalid substructure %r", substructure)
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
            detail=f"Invalid SMILES: {substructure!r}",
        )

    cache_key = f"subsearch:{substructure}"
    cached = get_cached_result(cache_key)
    if cached is not None:
        logger.debug("POST /search/substructure - Cache hit for %s", substructure)
        return {"task_id": None, "status": "completed", "result": cached}

    task = run_substructure_search.delay(substructure)
    logger.info("POST /search/substructure - Task %s queued", task.id)
    return {"task_id": task.id, "status": task.status}


@app.get("/search/substructure/{task_id}")
def get_substructure_task(task_id: str):
    logger.info("GET /search/substructure/%s - Checking task status", task_id)
    task_result = AsyncResult(task_id, app=celery_app)
    state = task_result.state

    if state == "PENDING":
        return {"task_id": task_id, "status": "pending"}
    if state == "STARTED":
        return {"task_id": task_id, "status": "in-progress"}
    if state == "FAILURE":
        logger.error("GET /search/substructure/%s - Task failed: %s", task_id, task_result.result)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(task_result.result),
        )
    if state == "SUCCESS":
        data = task_result.result or []
        return {"task_id": task_id, "status": "completed", "result": data}

    return {"task_id": task_id, "status": state.lower()}

