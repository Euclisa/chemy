from dataclasses import asdict, dataclass, field


@dataclass(frozen=True)
class ReactionTask:
    task_id: str
    cid: int
    compound_name: str
    preset: str
    position: str
    scope: str
    prompt: str
    existing_reactions: tuple = field(default_factory=tuple)
    complexity_band: str = None
    mode: str = "cold_start"
    metadata: dict = field(default_factory=dict)

    def to_dict(self):
        data = asdict(self)
        data["existing_reactions"] = list(self.existing_reactions)
        return data

    @classmethod
    def from_dict(cls, data):
        payload = data.copy()
        payload["existing_reactions"] = tuple(payload.get("existing_reactions") or [])
        return cls(**payload)


@dataclass(frozen=True)
class GenerationResult:
    task_id: str
    model: str
    raw_response: str
    reactions: tuple
    parse_stats: dict
    reviewed_raw_response: str = None
    reviewed_reactions: tuple = field(default_factory=tuple)
    reviewed_parse_stats: dict = field(default_factory=dict)
    usage: dict = field(default_factory=dict)
    latency_s: float = 0.0
    error: str = None

    def to_dict(self):
        data = asdict(self)
        data["reactions"] = list(self.reactions)
        data["reviewed_reactions"] = list(self.reviewed_reactions)
        return data


@dataclass(frozen=True)
class ReactionCandidate:
    candidate_id: str
    output_id: str
    task_id: str
    model: str
    reaction_index: int
    reaction: dict
    schema_valid: bool
    parsed: bool
    rid: str = None
    cid: int = None
    unmapped_names: tuple = field(default_factory=tuple)
    target_compound_present: bool = False
    target_position_correct: bool = False
    balance_success: bool = False
    duplicate_key: str = None

    def to_dict(self):
        data = asdict(self)
        data["unmapped_names"] = [list(item) for item in self.unmapped_names]
        return data


@dataclass(frozen=True)
class ValidationResult:
    sample_id: str
    reaction: dict
    validation_model: str
    validation_rounds: int
    validation_positives: int
    validation_confidence: float
    validation_votes: tuple
    rid: str = None
    cid: int = None
    task_id: str = None
    model: str = None

    def to_dict(self):
        data = asdict(self)
        data["validation_votes"] = list(self.validation_votes)
        return data


@dataclass(frozen=True)
class RepairResult:
    sample_id: str
    fixer_model: str
    original_reaction: dict
    fixed_reaction: dict = None
    repair_status: str = "null_repair"
    rid: str = None
    cid: int = None

    def to_dict(self):
        return asdict(self)
