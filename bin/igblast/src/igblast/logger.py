import io
import logging

LOGGER_FORMAT = '%(name)s | line %(lineno)-3d | %(levelname)-8s | %(message)s'


def set_logger(name: str, logger_format: str = LOGGER_FORMAT):
    """
    Initializes logging.

    :param name: name of the logging file
    :param logger_format: A logger format string
    """
    logger = logging.getLogger(name)
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter(logger_format))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    return logger


class TqdmToLogger(io.StringIO):
    """
    Output stream for TQDM which will output to logger module instead of
    the StdOut.
    """

    def __init__(self, logger, level=None):
        super().__init__()
        self.logger = logger
        self.buf = ''
        self.level = level or logging.INFO

    def write(self, buf: str) -> int:
        self.buf = buf.strip('\r\n\t ')
        return len(buf)

    def flush(self) -> None:
        self.logger.log(self.level, self.buf)
