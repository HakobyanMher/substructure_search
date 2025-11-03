import json
import redis
from os import getenv


redis_client = redis.from_url(
    getenv("REDIS_URL", "redis://redis:6379/0"),
    decode_responses=True,
)

def get_cached_result(key: str):
    data = redis_client.get(key)
    return json.loads(data) if data else None

def set_cache(key: str, value, expiration: int = 200):
    redis_client.setex(key, expiration, json.dumps(value))
