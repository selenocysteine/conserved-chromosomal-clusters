import logging
import sys

def init_logger():
    logger = logging.getLogger()
    if not logger.handlers:
        log_handler = logging.StreamHandler(sys.stdout)
        fmt = MyFormatter()
        log_handler.setFormatter(fmt)
        logger.addHandler(log_handler)

    return logger


class MyFormatter(logging.Formatter):
    # Ref:
    # https://stackoverflow.com/questions/14844970/modifying-logging-message-format-based-on-message-logging-level-in-python3
    standard_format = "%(asctime)s ## %(levelname)s: %(message)s"
    info_format = "%(asctime)s ## %(message)s"
    verbose_format = "%(asctime)s > %(message)s"

    def __init__(self):
        super().__init__(fmt="%(levelno)d: %(msg)s",
                         datefmt="%d-%m-%Y %H:%M:%S")

    def format(self, record):
        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.INFO:
            self._style._fmt = MyFormatter.info_format

        elif record.levelno == logging.VERBOSE:
            self._style._fmt = MyFormatter.verbose_format

        else:
            self._style._fmt = MyFormatter.standard_format

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result
