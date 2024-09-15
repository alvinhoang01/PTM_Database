from fastapi import FastAPI, UploadFile, File, Form
from pathlib import Path
import shutil
from datetime import datetime

app = FastAPI()

@app.post("/upload-fasta/")
async def upload_fasta(file: UploadFile = File(...), email: str = Form(...)):
    # Get the current year and month
    current_year = datetime.now().year
    current_month = datetime.now().strftime('%B')  # Gets month name like 'January', 'February'

    # Define the base directory
    base_dir = Path("C:\\Users\\maitr\\Documents\\nms\\PTM_Database_Copy\\ptmdatabase\\Storing_Fasta")

    # Create the folder structure: year/month/email
    user_dir = base_dir / str(current_year) / current_month / email
    user_dir.mkdir(parents=True, exist_ok=True)  # Create the directories if they don't exist

    # Define the path to save the file
    file_path = user_dir / file.filename

    # Save the uploaded file
    with file_path.open("wb") as buffer:
        shutil.copyfileobj(file.file, buffer)

    return {"filename": file.filename, "status": "uploaded", "saved_to": str(file_path)}
