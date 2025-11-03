from sqlalchemy import Column, String
from src.database import Base


class Molecule(Base):
    __tablename__ = "molecules"

    id = Column(String, primary_key=True, index=True)
    smiles = Column(String, nullable=False, index=True)

