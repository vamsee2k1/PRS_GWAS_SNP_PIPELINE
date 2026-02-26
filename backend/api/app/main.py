from __future__ import annotations

from fastapi import FastAPI

from app.deps import get_settings
from app.routes import health, jobs


def create_app() -> FastAPI:
    settings = get_settings()
    app = FastAPI(
        title=settings.app_name,
        version="0.1.0-mvp",
        description="Backend API for orchestrating pipeline jobs via Supabase + AWS VM.",
    )
    app.include_router(health.router)
    app.include_router(jobs.router)

    @app.get("/", tags=["root"])
    def root() -> dict[str, str]:
        return {"message": settings.app_name, "docs": "/docs"}

    return app


app = create_app()
