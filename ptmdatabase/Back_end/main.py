from fastapi import FastAPI, UploadFile, File, Form
from pathlib import Path
import shutil, re
from datetime import datetime

app = FastAPI()
VERSION = "upload-route-v3"

def slug(s: str) -> str:
    s = (s or "").strip().lower()
    s = re.sub(r'[^a-z0-9]+', '-', s)
    return s.strip('-') or "na"

@app.post("/upload-fasta/")
async def upload_fasta(
    file: UploadFile = File(...),
    username: str = Form(...),
    filename: str = Form(...),
    organization: str = Form(...),
    lab_pi: str = Form(...),
):
    now = datetime.now()
    year = f"{now:%Y}"
    month = now.strftime('%B')

    base_dir = Path(r"C:\Users\Administrator\Documents\Storing_Fasta")

    org_slug = slug(organization)
    pi_slug  = slug(lab_pi)
    user_slug = slug(username)

    # NEW tree
    user_dir = base_dir / year / month / f"org={org_slug}" / f"pi={pi_slug}" / f"user={user_slug}"
    user_dir.mkdir(parents=True, exist_ok=True)

    file_path = user_dir / filename

    print(f"[{VERSION}] org='{organization}' ({org_slug})  pi='{lab_pi}' ({pi_slug})  user='{username}' ({user_slug})", flush=True)
    print(f"[{VERSION}] SAVE TO: {file_path}", flush=True)

    with open(file_path, "wb") as out:
        shutil.copyfileobj(file.file, out)

    n_entries = sum(1 for line in open(file_path, "r", encoding="utf-8", errors="ignore") if line.startswith(">"))

    return {
        "status": "uploaded",
        "version": VERSION,
        "saved_to": str(file_path),
        "path_parts": {
            "base": str(base_dir),
            "year": year,
            "month": month,
            "org": f"org={org_slug}",
            "pi": f"pi={pi_slug}",
            "user": f"user={user_slug}",
            "filename": filename,
        },
        "echo": {"username": username, "organization": organization, "lab_pi": lab_pi},
        "n_entries": n_entries,
    }
