from fastapi import FastAPI, UploadFile, File, Form
from pathlib import Path
import shutil
from datetime import datetime

app = FastAPI()

@app.post("/upload-fasta/")
async def upload_fasta(file: UploadFile = File(...), username: str = Form(...), filename: str = Form(...)):
    # Get the current year and month
    current_year = datetime.now().year
    current_month = datetime.now().strftime('%B')

    # Define the base directory
    base_dir = Path("C:\\Users\\Administrator\\Documents\\Storing_Fasta")

    # Create the folder structure: year/month/username
    user_dir = base_dir / str(current_year) / current_month / username
    user_dir.mkdir(parents=True, exist_ok=True)  # Create the directories if they don't exist

    # Define the path to save the file, using the filename provided by the frontend
    file_path = user_dir / filename

    # Save the uploaded in-memory file to disk
    with open(file_path, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)

    return {"filename": filename, "status": "uploaded", "saved_to": str(file_path)}
