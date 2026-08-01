import os
from dotenv import load_dotenv

load_dotenv()

# Paths
RAW_DATA_PATH = os.getenv("RAW_DATA_PATH")
PROCESSED_DATA_PATH = os.getenv("PROCESSED_DATA_PATH")
OUTPUT_DIR = os.getenv("OUTPUT_DIR")

# Variables
TARGET = "hypertension"
DESIGN_COLS = ["cluster_id", "strata_id", "sample_weight"]

CATEGORICAL_VARS = [
    'age_group', 'sex', 'education_level', 'marital_status',
    'residence_type', 'division', 'wealth_index', 'BMI_category',
    'diabetes', 'household_size_cat', 'crowding_category',
    'electricity', 'improved_water', 'improved_toilet', 'mobile_phone'
]

# ML Config
RANDOM_SEED = int(os.getenv("RANDOM_SEED", 42))
TEST_SIZE = 0.2
VAL_SIZE = 0.25