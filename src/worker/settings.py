from lib.settings import Settings

CELERY_TASK_SERIALIZER = "pickle"
TASK_SERIALIZER = "pickle"
CELERY_ACCEPT_CONTENT = ["pickle"]
# CELERY_BROKER_URL = Settings.redis_uri
# CELERY_BACKEND_URL = Settings.redis_uri
BROKER_URL = Settings.redis_uri
BACKEND_URL = Settings.redis_uri
CELERYD_HIJACK_ROOT_LOGGER = False
