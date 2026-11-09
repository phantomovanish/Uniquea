# Uniquea in Public: GitHub, ORCID and What Comes Next

*Uniquea · November 9, 2026 · replaces the GitHub Release Plan of November 5 (kept for the record, House Rule 11)*

## Where things stand
| Item | Status |
|---|---|
| Repository | **github.com/phantomovanish/Uniquea**, public since November 5, 2026 |
| Author ID | **ORCID 0009-0005-8386-3999** (check digit verified) |
| Citation | CITATION.cff in the repository: GitHub shows a "Cite this repository" button crediting Stephen M. Griffin with the ORCID |
| Contents | Public Master (updated through November 9), EPR verdict, Tsirelson page, Hubble closure test, breathing limit, Orb checks, test 1 (SPARC), workbooks, figures, both audit packets |
| Safety | Every push is preceded by a private scan for withheld material (House Rule 16); the term list is never published |
| Owner's step | On the ORCID record, add github.com/phantomovanish/Uniquea as a website or work, so the link runs both ways |

## How these public sites work (Established)
- **GitHub** stores projects as *repositories* using **git**, a version-control system written by Linus Torvalds (2005). Every change is a *commit*: a dated snapshot with a note, and the whole history is kept, which is House Rule 11 in software. Anyone can read a public repository, copy it ("fork"), report problems ("issues") or suggest changes ("pull requests"); only the owner decides what goes in. GitHub was founded in 2008 and has been owned by Microsoft since 2018.
- **ORCID** (Open Researcher and Contributor ID) is a free, permanent ID for researchers, run by a nonprofit since 2012. It solves the "which Stephen Griffin?" problem. Journals, universities and archives attach it to work, and the owner controls what the record shows.
- **DOI** (Digital Object Identifier) is a permanent link to a specific version of a work. Papers and datasets get one so citations never break.
- **Zenodo**, run by CERN since 2013, gives any upload a DOI. It connects to GitHub: each GitHub "release" can be archived automatically with its own DOI, and it can sync to ORCID.
- **arXiv** is the main preprint server for physics and mathematics, founded by Paul Ginsparg in 1991 and now run by Cornell. Papers are moderated, not peer-reviewed. First-time authors need an **endorsement** from an established author in the category, and moderators may move a paper to a different category, such as general physics (physics.gen-ph).
- **LaTeX** is the typesetting language of mathematics and physics papers, Leslie Lamport's system (1984) built on Donald Knuth's TeX (1978). arXiv prefers LaTeX source.

## More public places, in order of value
1. **Zenodo** (recommended next): connect it to GitHub, make a release, and Uniquea gets a permanent DOI tied to the ORCID. No gatekeeper; the record is citable at once.
2. **Software Heritage**: archives public code repositories permanently; a GitHub repository can be saved there on request.
3. **OSF** (Open Science Framework, Center for Open Science) and **OSF Preprints**: projects, files and preprints with DOIs.
4. **arXiv**: the place physicists read, but it needs an endorsement and accepts moderation. The honest route is to submit a short, focused paper (for example, the SPARC test or the breathing limit) after an outside physicist has read it; the UCSD audit request is a natural path to an endorser.
5. **Google Scholar profile**: lists works once they are indexed (Zenodo or arXiv records count).
6. **ResearchGate**: a social site for researchers; useful for questions and contacts.
- *HAT note:* viXra accepts anything without moderation, and is not read by most physicists; it adds little credibility.

## A first taste of LaTeX
The One Math ladder, as it would appear in a paper:

```latex
\left(\frac{F}{F_\ast}\right)^{2^{n}} = \frac{m}{m_\ast}\left(\frac{a}{a_\ast}\right)^{2^{n+1}},
\qquad n = 0:\quad F = \frac{m\,a^{2}}{a_0}\;\;\text{(deep MOND; Milgrom 1983)}
```

and test 1's relation:

```latex
g_{\rm obs} = \sqrt{a_0\, g_{\rm bar}}, \qquad a_0 = (1.2\text{--}1.4)\times 10^{-10}\ {\rm m\,s^{-2}}
```

## Code, and Phantom's algorhymes
- **Code** has five meanings, and Uniquea uses all of them:
  1. a **cipher**: a hidden message (el-Bandito's poems, the SOS code 580/582);
  2. **computer code**: instructions a machine follows (every test script in the repository);
  3. the **genetic code**: DNA's three-letter words that spell 20 amino acids;
  4. a **legal code**: written rules (SMoG's department);
  5. a **moral code**: how we choose to live (HAT, the House Rules).
- **Algorithm** comes from the name of **Muhammad ibn Musa al-Khwarizmi** (9th century, Baghdad), whose book *al-jabr* also gave us the word **algebra** (House Rule 14).
- **Algorhyme**, Phantom's word, is algorithm + rhyme: a procedure with a beat. It has a famous ancestor: Radia Perlman's poem "Algorhyme" (1985), written to describe her spanning-tree algorithm, which still keeps computer networks from looping. Phantom turns everything into code because a code is a rule that can be checked: if the code runs, the idea holds; if it breaks, it says where.
