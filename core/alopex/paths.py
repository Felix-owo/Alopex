"""定位当前安装；不依赖工作目录或相邻 Pipeline。"""
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
CORE_ROOT = PIPELINE_ROOT / "core"
