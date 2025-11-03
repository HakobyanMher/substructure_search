from celery import Celery
from os import getenv


def _redis_url() -> str:
    return getenv("CELERY_BROKER_URL", getenv("REDIS_URL", "redis://redis:6379/0"))


celery_app = Celery(
    "molecule_store",
    broker=_redis_url(),
    backend=getenv("CELERY_RESULT_BACKEND", _redis_url()),
)

celery_app.conf.update(task_track_started=True)
celery_app.autodiscover_tasks(["src"])
