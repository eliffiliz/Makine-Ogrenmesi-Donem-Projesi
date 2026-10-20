"""Proje geneli sabitler. Tüm deneyler bu seed'i kullanır."""
from pathlib import Path

SEED = 42
SEEDS_FOR_REPEATED_CV = [42, 43, 44, 45, 46]  # final aşamasında en az 5 seed

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
GOLUB_DIR = DATA_DIR / "golub"
FIG_DIR = ROOT / "figures"

TRAIN_FILE = GOLUB_DIR / "data_set_ALL_AML_train.csv"
TEST_FILE = GOLUB_DIR / "data_set_ALL_AML_independent.csv"
LABEL_FILE = GOLUB_DIR / "actual.csv"

LABEL_MAP = {"ALL": 0, "AML": 1}
EXPECTED_TRAIN_COUNTS = {"ALL": 27, "AML": 11}
EXPECTED_TEST_COUNTS = {"ALL": 20, "AML": 14}
