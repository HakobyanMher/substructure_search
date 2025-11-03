from celery.utils.log import get_task_logger

from src.celery_app import celery_app
from src.cache import set_cache
from src.database import SessionLocal
from src.models import Molecule
from src.search import substructure_search

logger = get_task_logger(__name__)


@celery_app.task(name="tasks.run_substructure_search")
def run_substructure_search(substructure: str) -> list[dict[str, str]]:
    logger.info("Substructure search task started for %s", substructure)
    db = SessionLocal()
    try:
        molecules = db.query(Molecule).all()
        smiles_list = [m.smiles for m in molecules]
    finally:
        db.close()

    matches = set(substructure_search(smiles_list, substructure))
    result = [{"id": m.id, "smiles": m.smiles} for m in molecules if m.smiles in matches]
    set_cache(f"subsearch:{substructure}", result)
    logger.info("Substructure search task finished with %d results", len(result))
    return result
