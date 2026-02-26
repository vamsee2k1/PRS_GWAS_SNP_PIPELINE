from __future__ import annotations

from fastapi import APIRouter

from app.deps import get_settings

router = APIRouter(prefix="/health", tags=["health"])


@router.get("/live")
def live() -> dict[str, str]:
    return {"status": "ok"}


@router.get("/ready")
def ready() -> dict[str, str]:
    settings = get_settings()
    return {
        "status": "ok",
        "env": settings.app_env,
        "repository_backend": settings.repository_backend,
        "runner_backend": settings.runner_backend,
    }
