"""CLI에서 단계별로 보고할 수 있는 예외."""


class ResolverError(RuntimeError):
    error_code = "RESOLVER_ERROR"
    stage = "unknown"


class SequenceInputError(ResolverError):
    error_code = "INVALID_SEQUENCE"
    stage = "input"


class DatabaseUnavailableError(ResolverError):
    error_code = "DATABASE_UNAVAILABLE"
    stage = "database"


class StructureDownloadError(ResolverError):
    error_code = "STRUCTURE_DOWNLOAD_FAILED"
    stage = "download"


class PredictionUnavailableError(ResolverError):
    error_code = "PREDICTION_UNAVAILABLE"
    stage = "esmfold2"


class OutputExistsError(ResolverError):
    error_code = "OUTPUT_EXISTS"
    stage = "output"
