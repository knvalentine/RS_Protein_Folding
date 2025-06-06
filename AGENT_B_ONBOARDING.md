# Agent B Onboarding Guide

## Welcome Agent B - The Biophysicist!

### Quick Start
1. **Repository Location**: `/Users/kelsey/Documents/RS_Protein_Folding`
2. **GitHub**: https://github.com/knvalentine/RS_Protein_Folding

### How to Communicate with Agent A

#### Option 1: Direct File Updates (Recommended)
Edit these files directly in the repository:
- `AGENT_COMMUNICATIONS.md` - Main communication hub
- `AGENT_B_FINDINGS.md` - Your dedicated findings document (create as needed)

#### Option 2: Through Code Comments
Add comments in relevant code files with the tag `# AGENT_B:` for easy searching

#### Option 3: Create New Analysis Files
Create new files in `src/analysis/` for your biophysical analyses

### First Steps for Agent B

1. **Confirm Communication**:
   - Edit `AGENT_COMMUNICATIONS.md`
   - Update the "Agent B Status" section with current timestamp
   - Add a test message in the findings section

2. **Review Key Files**:
   ```
   docs/ROADMAP_CLEAN.md          # Project history
   CURRENT_APPROACH_SUMMARY.md    # Where we're stuck
   src/tests/test_three_proteins_quick.py  # See the problem
   ```

3. **Run a Test**:
   ```bash
   cd /Users/kelsey/Documents/RS_Protein_Folding
   cd src/tests
   python3 test_three_proteins_quick.py
   ```

### Communication Test Checklist

- [ ] Can access the repository at the path above
- [ ] Can edit AGENT_COMMUNICATIONS.md
- [ ] Can run the test script
- [ ] Understands the β-sheet prediction problem

### Key Questions to Start

1. Look at the test output - why are Trp-cage and BBA5 folding so fast?
2. What's different about Villin (which works well)?
3. What template features could predict β-sheet formation?

### Your Mission
Help us predict R (β-registries) from template-stage information, so our folding time predictions match experiments.

---
*Agent A is monitoring this file and AGENT_COMMUNICATIONS.md for updates* 