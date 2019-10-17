def return_config(verbose):
    logging_config = {"version": 1, "disable_existing_loggers": False}
    fmt = '[%(levelname)-5.5s] %(asctime)s :: => %(message)s'
    if verbose:
        fmt = '[%(levelname)-5.5s] %(asctime)s %(pathname)-80.80s :: %(lineno)s => %(message)s'
    logging_config["formatters"] = {"basic": {"format": fmt, "datefmt": "%Y-%m-%d %H:%M:%S"}}

    logging_config["handlers"] = {
            "console": {
                "class": "logging.StreamHandler",
                "level": "INFO",
                "formatter": "basic",
                "stream": "ext://sys.stdout"
            },
            "info_file_handler": {
                "class": "logging.handlers.RotatingFileHandler",
                "level": "INFO",
                "formatter": "basic",
                "filename": "log.txt",
                "maxBytes": 10485760,
                "backupCount": 20,
                "encoding": "utf8"
            }
        }
    libLevel = "WARNING"
    if verbose:
        libLevel = "INFO"
    logging_config["loggers"] = {
            "werkzeug": {
                "level": libLevel,
                "handlers": ["console", "info_file_handler"],
                "propagate": False
            }
        }
    logging_config["root"] = {
            "level": "INFO",
            "handlers": ["console", "info_file_handler"]
        }

    return logging_config


