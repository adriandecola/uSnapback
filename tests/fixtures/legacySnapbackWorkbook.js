// Frozen values transcribed from Carl's local reference workbooks.
// Keep source locations and limitations with the data so these historical
// outputs are never mistaken for fully specified scientific ground truth.
export const LEGACY_WORKBOOK_SOURCE = Object.freeze({
  "input": "../Snapback Tm for Adrian/Autocalc Datafile.xlsx :: Autocalc Datafile to 20!A1:I48",
  "inputSha256": "37da7ae010e29118666534c6f2b2e3aaccb4e8acc13c104c706b99b71f6c8ff5",
  "output": "../Snapback Tm for Adrian/Output Data - Adrian.xlsx :: Output 20!A1:K48",
  "outputSha256": "c585761777e2e3cce21c3f8f9c590b2ac953988dc28092464c68d038e96bf031",
  "hughRefit": "../Snapback Tm for Adrian/Output Data - Adrian with Hugh Empirical Refit Formula.xlsx :: Hugh Refit!A1:S45",
  "hughRefitSha256": "530d9378b99d56a3fa1ff804ac99e838783fcb2c058200e2e9870f40b39ece9d",
  "limitations": "The files omit effective monovalent salt, strand concentrations, free-vs-total Mg, complete terminal tetrads, and loop-end construction details."
});

export const UNRECONSTRUCTABLE_END_ROWS = Object.freeze({
  8: "The recorded SNP3 left end is Watson-Crick complementary, so the identity of the engineered loop-side mismatch is missing.",
  9: "The recorded SNP3 left end is Watson-Crick complementary, so the identity of the engineered loop-side mismatch is missing."
});

export const LEGACY_DELTA_TM_SUMMARY = Object.freeze({
  empiricalNoEnds: Object.freeze({ n: 46, bias: 0.08695652173912982, mae: 1.617391304347828, rmse: 2.439885955994843, sampleSd: 2.465279672721596 }),
  empiricalWithEnds: Object.freeze({ n: 46, bias: -0.3739130434782614, mae: 1.5347826086956515, rmse: 2.0186090777304173, sampleSd: 2.005596035830045 }),
  santaLucia: Object.freeze({ n: 46, bias: -0.6217391304347827, mae: 1.665217391304347, rmse: 2.0913902178881534, sampleSd: 2.0189010259685842 }),
  rochester: Object.freeze({ n: 46, bias: -0.3434782608695652, mae: 1.61304347826087, rmse: 1.9554828159506434, sampleSd: 1.9463529626130152 })
});

export const LEGACY_ALL_PEAK_SUMMARY = Object.freeze({
  empiricalNoEnds: Object.freeze({ n: 93, bias: 0.43978494623656017, mae: 3.321505376344088, sampleSd: 4.361121210665756, rSquared: 0.7036 }),
  empiricalWithEnds: Object.freeze({ n: 93, bias: 3.091397849462367, mae: 4.121505376344088, sampleSd: 3.789535244934928, rSquared: 0.7383 }),
  santaLucia: Object.freeze({ n: 93, bias: -1.4881720430107515, mae: 2.776344086021506, sampleSd: 3.2728411225658, rSquared: 0.7817 }),
  rochester: Object.freeze({ n: 93, bias: -10.273118279569891, mae: 10.273118279569891, sampleSd: 4.001063439749668, rSquared: 0.722 })
});

export const LEGACY_WORKBOOK_ROWS = Object.freeze([
  {
    "workbookRow": 2,
    "snpName": "SNP1",
    "leftMismatch": "G",
    "coreSequence": "GCCAGCCAGTT(CT)TGCAAAGAACC",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 4,
    "loopLength": 51,
    "operator": "Hugh",
    "measured": {
      "matched": 77,
      "mismatched": 72.3
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 79.4,
        "mismatched": 75.3
      },
      "empiricalWithEnds": {
        "matched": 80.4,
        "mismatched": 76.6
      },
      "santaLucia": {
        "matched": 75.9,
        "mismatched": 71.9
      },
      "rochester": {
        "matched": 69.5,
        "mismatched": 65.5
      }
    }
  },
  {
    "workbookRow": 3,
    "snpName": "SNP1",
    "leftMismatch": "G",
    "coreSequence": "GCCAGCCAGTT(CT)TGCAAAGAACC",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 2,
    "loopLength": 51,
    "operator": "Ross",
    "measured": {
      "matched": 73.6,
      "mismatched": 69.1
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 77.8,
        "mismatched": 73.6
      },
      "empiricalWithEnds": {
        "matched": 79,
        "mismatched": 75.1
      },
      "santaLucia": {
        "matched": 74,
        "mismatched": 70
      },
      "rochester": {
        "matched": 67.9,
        "mismatched": 63.8
      }
    }
  },
  {
    "workbookRow": 4,
    "snpName": "SNP1-7C",
    "leftMismatch": "C",
    "coreSequence": "CCAGCCAGTT(CT)TGCAAAG",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 4,
    "loopLength": 53,
    "operator": "Hugh",
    "measured": {
      "matched": 73.7,
      "mismatched": 67.2
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 73.5,
        "mismatched": 67.2
      },
      "empiricalWithEnds": {
        "matched": 76,
        "mismatched": 69.9
      },
      "santaLucia": {
        "matched": 73.7,
        "mismatched": 67.6
      },
      "rochester": {
        "matched": 65,
        "mismatched": 58.9
      }
    }
  },
  {
    "workbookRow": 5,
    "snpName": "SNP1",
    "leftMismatch": "T",
    "coreSequence": "TTCTTTGCA(GA)AACTGGCTGG",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3,
    "loopLength": 41,
    "operator": "Adrian",
    "measured": {
      "matched": 74.2,
      "mismatched": 67
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 75.9,
        "mismatched": 69.7
      },
      "empiricalWithEnds": {
        "matched": 78.9,
        "mismatched": 73.3
      },
      "santaLucia": {
        "matched": 75.7,
        "mismatched": 71.1
      },
      "rochester": {
        "matched": 69.7,
        "mismatched": 64.6
      }
    }
  },
  {
    "workbookRow": 6,
    "snpName": "SNP2",
    "leftMismatch": "T",
    "coreSequence": "GACAAAAAAT(GA)TACTTATACTT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 3,
    "loopLength": 24,
    "operator": "Hugh",
    "measured": {
      "matched": 68.2,
      "mismatched": 60.9
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 67.7,
        "mismatched": 60.9
      },
      "empiricalWithEnds": {
        "matched": 69.5,
        "mismatched": 63.1
      },
      "santaLucia": {
        "matched": 64.4,
        "mismatched": 58.8
      },
      "rochester": {
        "matched": 61.9,
        "mismatched": 56.5
      }
    }
  },
  {
    "workbookRow": 7,
    "snpName": "SNP2",
    "leftMismatch": "T",
    "coreSequence": "GACAAAAAAT(GA)TACTTATACTT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 2,
    "loopLength": 24,
    "operator": "Ross",
    "measured": {
      "matched": 69.4,
      "mismatched": 61.7
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 66.4,
        "mismatched": 59.7
      },
      "empiricalWithEnds": {
        "matched": 68.3,
        "mismatched": 61.9
      },
      "santaLucia": {
        "matched": 62.9,
        "mismatched": 57.3
      },
      "rochester": {
        "matched": 60.5,
        "mismatched": 55.2
      }
    }
  },
  {
    "workbookRow": 8,
    "snpName": "SNP3",
    "leftMismatch": "T",
    "coreSequence": "ACATGACCAAGC(TA)CCCCAAATTCA",
    "rightMismatch": "A",
    "reportedMagnesiumMm": 3.5,
    "loopLength": 88,
    "operator": "Hugh",
    "measured": {
      "matched": 73.6,
      "mismatched": 69.1
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 77.7,
        "mismatched": 74.1
      },
      "empiricalWithEnds": {
        "matched": 80,
        "mismatched": 76.8
      },
      "santaLucia": {
        "matched": 73.9,
        "mismatched": 70.4
      },
      "rochester": {
        "matched": 62.7,
        "mismatched": 58.9
      }
    }
  },
  {
    "workbookRow": 9,
    "snpName": "SNP3",
    "leftMismatch": "T",
    "coreSequence": "ACATGACCAAGC(TA)CCCCAAATTCA",
    "rightMismatch": "A",
    "reportedMagnesiumMm": 2,
    "loopLength": 88,
    "operator": "Ross",
    "measured": {
      "matched": 74,
      "mismatched": 69.3
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 76.4,
        "mismatched": 72.8
      },
      "empiricalWithEnds": {
        "matched": 78.7,
        "mismatched": 75.5
      },
      "santaLucia": {
        "matched": 72.3,
        "mismatched": 68.8
      },
      "rochester": {
        "matched": 61.2,
        "mismatched": 57.5
      }
    }
  },
  {
    "workbookRow": 10,
    "snpName": "SNP4",
    "leftMismatch": "T",
    "coreSequence": "CACCTTCTA(TA)AATTTTAA",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3.5,
    "loopLength": 38,
    "operator": "Hugh",
    "measured": {
      "matched": 60.3,
      "mismatched": 53.6
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 60,
        "mismatched": 54.2
      },
      "empiricalWithEnds": {
        "matched": 62.6,
        "mismatched": 57.1
      },
      "santaLucia": {
        "matched": 58.8,
        "mismatched": 54.8
      },
      "rochester": {
        "matched": 54.2,
        "mismatched": 50
      }
    }
  },
  {
    "workbookRow": 11,
    "snpName": "SNP4",
    "leftMismatch": "T",
    "coreSequence": "CACCTTCTA(TA)AATTTTAA",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3,
    "loopLength": 38,
    "operator": "Ross",
    "measured": {
      "matched": 62.3,
      "mismatched": 57.1
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 59.5,
        "mismatched": 53.7
      },
      "empiricalWithEnds": {
        "matched": 62.1,
        "mismatched": 56.6
      },
      "santaLucia": {
        "matched": 58.2,
        "mismatched": 54.2
      },
      "rochester": {
        "matched": 53.8,
        "mismatched": 49.6
      }
    }
  },
  {
    "workbookRow": 12,
    "snpName": "SNP5",
    "leftMismatch": "A",
    "coreSequence": "ACGTTCAC(CT)ACACAGTTA",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3.5,
    "loopLength": 90,
    "operator": "Hugh",
    "measured": {
      "matched": 68.2,
      "mismatched": 63.1
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 71.5,
        "mismatched": 65.6
      },
      "empiricalWithEnds": {
        "matched": 73.8,
        "mismatched": 68.3
      },
      "santaLucia": {
        "matched": 68.9,
        "mismatched": 63.1
      },
      "rochester": {
        "matched": 54.3,
        "mismatched": 48.2
      }
    }
  },
  {
    "workbookRow": 13,
    "snpName": "SNP5",
    "leftMismatch": "A",
    "coreSequence": "ACGTTCAC(CT)ACACAGTTA",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3,
    "loopLength": 90,
    "operator": "Ross",
    "measured": {
      "matched": 70.6,
      "mismatched": 64.3
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 71.1,
        "mismatched": 65.2
      },
      "empiricalWithEnds": {
        "matched": 73.4,
        "mismatched": 67.9
      },
      "santaLucia": {
        "matched": 68.4,
        "mismatched": 62.6
      },
      "rochester": {
        "matched": 53.9,
        "mismatched": 47.8
      }
    }
  },
  {
    "workbookRow": 14,
    "snpName": "SNP6",
    "leftMismatch": "C",
    "coreSequence": "TTTTTTTT(AT)AAATACCATCT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 3.5,
    "loopLength": 79,
    "operator": "Hugh",
    "measured": {
      "matched": 59.4,
      "mismatched": 55.8
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 61.5,
        "mismatched": 56.5
      },
      "empiricalWithEnds": {
        "matched": 63.5,
        "mismatched": 58.7
      },
      "santaLucia": {
        "matched": 56.5,
        "mismatched": 51.7
      },
      "rochester": {
        "matched": 46,
        "mismatched": 40.7
      }
    }
  },
  {
    "workbookRow": 15,
    "snpName": "SNP6",
    "leftMismatch": "C",
    "coreSequence": "TTTTTTTT(AT)AAATACCATCT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 4,
    "loopLength": 79,
    "operator": "Hugh",
    "measured": {
      "matched": 59.9,
      "mismatched": 56.4
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 61.9,
        "mismatched": 56.8
      },
      "empiricalWithEnds": {
        "matched": 63.8,
        "mismatched": 59
      },
      "santaLucia": {
        "matched": 56.9,
        "mismatched": 52
      },
      "rochester": {
        "matched": 46.3,
        "mismatched": 41.1
      }
    }
  },
  {
    "workbookRow": 16,
    "snpName": "SNP6",
    "leftMismatch": "C",
    "coreSequence": "TTTTTTTT(AT)AAATACCATCT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 3,
    "loopLength": 79,
    "operator": "Ross",
    "measured": {
      "matched": 61.1,
      "mismatched": 57.9
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 61.1,
        "mismatched": 56.1
      },
      "empiricalWithEnds": {
        "matched": 63.1,
        "mismatched": 58.3
      },
      "santaLucia": {
        "matched": 56,
        "mismatched": 51.2
      },
      "rochester": {
        "matched": 45.6,
        "mismatched": 40.3
      }
    }
  },
  {
    "workbookRow": 17,
    "snpName": "SNP7",
    "leftMismatch": "C",
    "coreSequence": "CAACGAGC(GA)TCTTGTAAAAT",
    "rightMismatch": "T",
    "reportedMagnesiumMm": 3,
    "loopLength": 42,
    "operator": "Hugh",
    "measured": {
      "matched": 71.2,
      "mismatched": 61.4
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 72.9,
        "mismatched": 65.4
      },
      "empiricalWithEnds": {
        "matched": 75.1,
        "mismatched": 67.9
      },
      "santaLucia": {
        "matched": 71.4,
        "mismatched": 64.9
      },
      "rochester": {
        "matched": 65.6,
        "mismatched": 58.8
      }
    }
  },
  {
    "workbookRow": 18,
    "snpName": "SNP7",
    "leftMismatch": "C",
    "coreSequence": "CAACGAGC(GA)TCTTGTAAAAT",
    "rightMismatch": "T",
    "reportedMagnesiumMm": 2,
    "loopLength": 42,
    "operator": "Ross",
    "measured": {
      "matched": 73.6,
      "mismatched": 64.1
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 71.6,
        "mismatched": 64.2
      },
      "empiricalWithEnds": {
        "matched": 73.9,
        "mismatched": 66.7
      },
      "santaLucia": {
        "matched": 69.9,
        "mismatched": 63.4
      },
      "rochester": {
        "matched": 64.3,
        "mismatched": 57.6
      }
    }
  },
  {
    "workbookRow": 19,
    "snpName": "SNP7",
    "leftMismatch": "A",
    "coreSequence": "CGAGC(GA)TCTTGT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 3,
    "loopLength": 44,
    "operator": "Adrian",
    "measured": {
      "matched": 65.9,
      "mismatched": 56.3
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 58.3,
        "mismatched": 37.5
      },
      "empiricalWithEnds": {
        "matched": 61.8,
        "mismatched": 44.5
      },
      "santaLucia": {
        "matched": 61.2,
        "mismatched": 45.4
      },
      "rochester": {
        "matched": 53.9,
        "mismatched": 39.2
      }
    }
  },
  {
    "workbookRow": 20,
    "snpName": "SNP8",
    "leftMismatch": "A",
    "coreSequence": "AATGTAATGCA(GA)TTGAGTTAGAAAT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 3.5,
    "loopLength": 74,
    "operator": "Hugh",
    "measured": {
      "matched": 68.2,
      "mismatched": 60.4
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 72.8,
        "mismatched": 68.4
      },
      "empiricalWithEnds": {
        "matched": 75.5,
        "mismatched": 71.2
      },
      "santaLucia": {
        "matched": 69.3,
        "mismatched": 65.3
      },
      "rochester": {
        "matched": 60.1,
        "mismatched": 55.3
      }
    }
  },
  {
    "workbookRow": 21,
    "snpName": "SNP8",
    "leftMismatch": "A",
    "coreSequence": "AATGTAATGCA(GA)TTGAGTTAGAAAT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 2,
    "loopLength": 74,
    "operator": "Ross",
    "measured": {
      "matched": 69,
      "mismatched": 62
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 71.4,
        "mismatched": 66.9
      },
      "empiricalWithEnds": {
        "matched": 74,
        "mismatched": 69.8
      },
      "santaLucia": {
        "matched": 67.6,
        "mismatched": 63.5
      },
      "rochester": {
        "matched": 58.5,
        "mismatched": 53.7
      }
    }
  },
  {
    "workbookRow": 22,
    "snpName": "SNP8",
    "leftMismatch": "A",
    "coreSequence": "AATGCA(GA)TTGAGTT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 4,
    "loopLength": 79,
    "operator": "Adrian",
    "measured": {
      "matched": null,
      "mismatched": 48
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 61.5,
        "mismatched": 47.9
      },
      "empiricalWithEnds": {
        "matched": 64.7,
        "mismatched": 51.7
      },
      "santaLucia": {
        "matched": 62.3,
        "mismatched": 50.6
      },
      "rochester": {
        "matched": 46.1,
        "mismatched": 33
      }
    }
  },
  {
    "workbookRow": 23,
    "snpName": "SNP9",
    "leftMismatch": "A",
    "coreSequence": "CACTGGGACAG(CG)AATAAAACACCT",
    "rightMismatch": "T",
    "reportedMagnesiumMm": 3.5,
    "loopLength": 71,
    "operator": "Hugh",
    "measured": {
      "matched": 73,
      "mismatched": 66.6
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 77.5,
        "mismatched": 74.3
      },
      "empiricalWithEnds": {
        "matched": 79.2,
        "mismatched": 76.2
      },
      "santaLucia": {
        "matched": 73.7,
        "mismatched": 70.7
      },
      "rochester": {
        "matched": 64.7,
        "mismatched": 61.3
      }
    }
  },
  {
    "workbookRow": 24,
    "snpName": "SNP9",
    "leftMismatch": "A",
    "coreSequence": "CACTGGGACAG(CG)AATAAAACACCT",
    "rightMismatch": "T",
    "reportedMagnesiumMm": 3,
    "loopLength": 71,
    "operator": "Ross",
    "measured": {
      "matched": 74.6,
      "mismatched": 69.2
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 77.2,
        "mismatched": 74
      },
      "empiricalWithEnds": {
        "matched": 78.9,
        "mismatched": 75.9
      },
      "santaLucia": {
        "matched": 73.3,
        "mismatched": 70.3
      },
      "rochester": {
        "matched": 64.3,
        "mismatched": 60.9
      }
    }
  },
  {
    "workbookRow": 25,
    "snpName": "SNP9",
    "leftMismatch": "T",
    "coreSequence": "TTTATT(CG)CTGTCCC",
    "rightMismatch": "T",
    "reportedMagnesiumMm": 3,
    "loopLength": 34,
    "operator": "Adrian",
    "measured": {
      "matched": 63.3,
      "mismatched": 52.9
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 58.6,
        "mismatched": 51.1
      },
      "empiricalWithEnds": {
        "matched": 62.1,
        "mismatched": 54.8
      },
      "santaLucia": {
        "matched": 61.4,
        "mismatched": 55
      },
      "rochester": {
        "matched": 56.4,
        "mismatched": 50.5
      }
    }
  },
  {
    "workbookRow": 26,
    "snpName": "SNP10",
    "leftMismatch": "G",
    "coreSequence": "GGGAAG(TC)CAGCAGG",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3,
    "loopLength": 50,
    "operator": "Hugh",
    "measured": {
      "matched": 67.5,
      "mismatched": 61.8
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 63.7,
        "mismatched": 52.8
      },
      "empiricalWithEnds": {
        "matched": 67.7,
        "mismatched": 58.8
      },
      "santaLucia": {
        "matched": 65.8,
        "mismatched": 57.8
      },
      "rochester": {
        "matched": 57.5,
        "mismatched": 49.3
      }
    }
  },
  {
    "workbookRow": 27,
    "snpName": "SNP10",
    "leftMismatch": "G",
    "coreSequence": "GGGAAG(TC)CAGCAGG",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 2,
    "loopLength": 50,
    "operator": "Ross",
    "measured": {
      "matched": 70.5,
      "mismatched": 64
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 62.1,
        "mismatched": 51.3
      },
      "empiricalWithEnds": {
        "matched": 66.3,
        "mismatched": 57.5
      },
      "santaLucia": {
        "matched": 64,
        "mismatched": 56.2
      },
      "rochester": {
        "matched": 56,
        "mismatched": 47.9
      }
    }
  },
  {
    "workbookRow": 28,
    "snpName": "SNP10",
    "leftMismatch": "C",
    "coreSequence": "CCTGCTG(TC)CTTCCCGT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 4,
    "loopLength": 46,
    "operator": "Adrian",
    "measured": {
      "matched": 77.9,
      "mismatched": 68.5
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 72.7,
        "mismatched": 65.2
      },
      "empiricalWithEnds": {
        "matched": 74.3,
        "mismatched": 67.6
      },
      "santaLucia": {
        "matched": 72.6,
        "mismatched": 66.8
      },
      "rochester": {
        "matched": 65,
        "mismatched": 58.8
      }
    }
  },
  {
    "workbookRow": 29,
    "snpName": "SNP11",
    "leftMismatch": "G",
    "coreSequence": "TAGCTATCAAC(GA)GGGCTGATATTCA",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3,
    "loopLength": 73,
    "operator": "Hugh",
    "measured": {
      "matched": 73.6,
      "mismatched": 68.3
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 78,
        "mismatched": 71.6
      },
      "empiricalWithEnds": {
        "matched": 78.4,
        "mismatched": 72.4
      },
      "santaLucia": {
        "matched": 72.2,
        "mismatched": 65.7
      },
      "rochester": {
        "matched": 63.7,
        "mismatched": 57
      }
    }
  },
  {
    "workbookRow": 30,
    "snpName": "SNP11",
    "leftMismatch": "G",
    "coreSequence": "TAGCTATCAAC(GA)GGGCTGATATTCA",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3,
    "loopLength": 73,
    "operator": "Ross",
    "measured": {
      "matched": 76.5,
      "mismatched": 70.2
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 78,
        "mismatched": 71.6
      },
      "empiricalWithEnds": {
        "matched": 78.4,
        "mismatched": 72.4
      },
      "santaLucia": {
        "matched": 72.2,
        "mismatched": 65.7
      },
      "rochester": {
        "matched": 63.7,
        "mismatched": 57
      }
    }
  },
  {
    "workbookRow": 31,
    "snpName": "SNP12",
    "leftMismatch": "C",
    "coreSequence": "TTGGAGAAAAG(TC)AGACAATATGATTT",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3,
    "loopLength": 83,
    "operator": "Hugh",
    "measured": {
      "matched": 65.3,
      "mismatched": 61.3
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 71.7,
        "mismatched": 68.5
      },
      "empiricalWithEnds": {
        "matched": 73.9,
        "mismatched": 70.8
      },
      "santaLucia": {
        "matched": 67.1,
        "mismatched": 64
      },
      "rochester": {
        "matched": 57.1,
        "mismatched": 53.4
      }
    }
  },
  {
    "workbookRow": 32,
    "snpName": "SNP12",
    "leftMismatch": "C",
    "coreSequence": "TTGGAGAAAAG(TC)AGACAATATGATTT",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 2,
    "loopLength": 83,
    "operator": "Ross",
    "measured": {
      "matched": 68.4,
      "mismatched": 63.7
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 70.6,
        "mismatched": 67.4
      },
      "empiricalWithEnds": {
        "matched": 72.8,
        "mismatched": 69.8
      },
      "santaLucia": {
        "matched": 65.8,
        "mismatched": 62.7
      },
      "rochester": {
        "matched": 56,
        "mismatched": 52.3
      }
    }
  },
  {
    "workbookRow": 33,
    "snpName": "SNP13",
    "leftMismatch": "T",
    "coreSequence": "CAACACTTGC(GA)AATATTGCT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 3,
    "loopLength": 91,
    "operator": "Hugh",
    "measured": {
      "matched": 66.3,
      "mismatched": 58.4
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 70.4,
        "mismatched": 62.2
      },
      "empiricalWithEnds": {
        "matched": 71.6,
        "mismatched": 63.9
      },
      "santaLucia": {
        "matched": 65.9,
        "mismatched": 58.1
      },
      "rochester": {
        "matched": 52.7,
        "mismatched": 43.5
      }
    }
  },
  {
    "workbookRow": 34,
    "snpName": "SNP13",
    "leftMismatch": "T",
    "coreSequence": "CAACACTTGC(GA)AATATTGCT",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 3,
    "loopLength": 91,
    "operator": "Ross",
    "measured": {
      "matched": 69.2,
      "mismatched": 61.9
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 70.4,
        "mismatched": 62.2
      },
      "empiricalWithEnds": {
        "matched": 71.6,
        "mismatched": 63.9
      },
      "santaLucia": {
        "matched": 65.9,
        "mismatched": 58.1
      },
      "rochester": {
        "matched": 52.7,
        "mismatched": 43.5
      }
    }
  },
  {
    "workbookRow": 35,
    "snpName": "SNP14",
    "leftMismatch": "T",
    "coreSequence": "GAAACAACA(TA)GCAGCTCAGGGTG",
    "rightMismatch": "A",
    "reportedMagnesiumMm": 3,
    "loopLength": 74,
    "operator": "Hugh",
    "measured": {
      "matched": 74.2,
      "mismatched": 69.5
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 78.8,
        "mismatched": 75.3
      },
      "empiricalWithEnds": {
        "matched": 81.5,
        "mismatched": 78.2
      },
      "santaLucia": {
        "matched": 77,
        "mismatched": 73.8
      },
      "rochester": {
        "matched": 66.4,
        "mismatched": 62.6
      }
    }
  },
  {
    "workbookRow": 36,
    "snpName": "SNP14",
    "leftMismatch": "T",
    "coreSequence": "GAAACAACA(TA)GCAGCTCAGGGTG",
    "rightMismatch": "A",
    "reportedMagnesiumMm": 2,
    "loopLength": 74,
    "operator": "Ross",
    "measured": {
      "matched": 75.4,
      "mismatched": 71.1
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 77.6,
        "mismatched": 74.2
      },
      "empiricalWithEnds": {
        "matched": 80.4,
        "mismatched": 77.2
      },
      "santaLucia": {
        "matched": 75.7,
        "mismatched": 72.5
      },
      "rochester": {
        "matched": 65.3,
        "mismatched": 61.6
      }
    }
  },
  {
    "workbookRow": 37,
    "snpName": "SNP15",
    "leftMismatch": "G",
    "coreSequence": "ATGTGC(TC)AGCTGAATGT",
    "rightMismatch": "T",
    "reportedMagnesiumMm": 3,
    "loopLength": 37,
    "operator": "Hugh",
    "measured": {
      "matched": 71.2,
      "mismatched": 65.5
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 69.4,
        "mismatched": 63
      },
      "empiricalWithEnds": {
        "matched": 74.3,
        "mismatched": 68.3
      },
      "santaLucia": {
        "matched": 72.4,
        "mismatched": 67.8
      },
      "rochester": {
        "matched": 66.3,
        "mismatched": 61.4
      }
    }
  },
  {
    "workbookRow": 38,
    "snpName": "SNP15",
    "leftMismatch": "G",
    "coreSequence": "ATGTGC(TC)AGCTGAATGT",
    "rightMismatch": "T",
    "reportedMagnesiumMm": 2,
    "loopLength": 37,
    "operator": "Ross",
    "measured": {
      "matched": 73.1,
      "mismatched": 68.6
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 67.9,
        "mismatched": 61.7
      },
      "empiricalWithEnds": {
        "matched": 72.9,
        "mismatched": 67
      },
      "santaLucia": {
        "matched": 70.7,
        "mismatched": 66.1
      },
      "rochester": {
        "matched": 64.8,
        "mismatched": 60
      }
    }
  },
  {
    "workbookRow": 39,
    "snpName": "SNP16",
    "leftMismatch": "C",
    "coreSequence": "ACATACATTC(GT)GTCTGCAGAG",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3,
    "loopLength": 82,
    "operator": "Hugh",
    "measured": {
      "matched": 68.6,
      "mismatched": 60.3
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 73.6,
        "mismatched": 65.8
      },
      "empiricalWithEnds": {
        "matched": 76.2,
        "mismatched": 68.6
      },
      "santaLucia": {
        "matched": 71,
        "mismatched": 62.6
      },
      "rochester": {
        "matched": 59.1,
        "mismatched": 50.8
      }
    }
  },
  {
    "workbookRow": 40,
    "snpName": "SNP16",
    "leftMismatch": "C",
    "coreSequence": "ACATACATTC(GT)GTCTGCAGAG",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3,
    "loopLength": 82,
    "operator": "Ross",
    "measured": {
      "matched": 72.3,
      "mismatched": 63.5
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 73.6,
        "mismatched": 65.8
      },
      "empiricalWithEnds": {
        "matched": 76.2,
        "mismatched": 68.6
      },
      "santaLucia": {
        "matched": 71,
        "mismatched": 62.6
      },
      "rochester": {
        "matched": 59.1,
        "mismatched": 50.8
      }
    }
  },
  {
    "workbookRow": 41,
    "snpName": "SNP17",
    "leftMismatch": "C",
    "coreSequence": "TTTTTACAAA(TC)ATATGATTATA",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 3,
    "loopLength": 46,
    "operator": "Hugh",
    "measured": {
      "matched": 64.2,
      "mismatched": 59.3
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 61.6,
        "mismatched": 57.2
      },
      "empiricalWithEnds": {
        "matched": 63.7,
        "mismatched": 59.5
      },
      "santaLucia": {
        "matched": 56.8,
        "mismatched": 53.5
      },
      "rochester": {
        "matched": 52.2,
        "mismatched": 48.6
      }
    }
  },
  {
    "workbookRow": 42,
    "snpName": "SNP17",
    "leftMismatch": "C",
    "coreSequence": "TTTTTACAAA(TC)ATATGATTATA",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 3,
    "loopLength": 46,
    "operator": "Ross",
    "measured": {
      "matched": 66.9,
      "mismatched": 62.3
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 61.6,
        "mismatched": 57.2
      },
      "empiricalWithEnds": {
        "matched": 63.7,
        "mismatched": 59.5
      },
      "santaLucia": {
        "matched": 56.8,
        "mismatched": 53.5
      },
      "rochester": {
        "matched": 52.2,
        "mismatched": 48.6
      }
    }
  },
  {
    "workbookRow": 43,
    "snpName": "SNP18",
    "leftMismatch": "T",
    "coreSequence": "TGTTTATTTTG(GT)GGCCTAGC",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 3,
    "loopLength": 39,
    "operator": "Hugh",
    "measured": {
      "matched": 71.3,
      "mismatched": 65.6
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 71.8,
        "mismatched": 62.8
      },
      "empiricalWithEnds": {
        "matched": 73.8,
        "mismatched": 65.4
      },
      "santaLucia": {
        "matched": 69.3,
        "mismatched": 60.6
      },
      "rochester": {
        "matched": 64.5,
        "mismatched": 56.2
      }
    }
  },
  {
    "workbookRow": 44,
    "snpName": "SNP18",
    "leftMismatch": "T",
    "coreSequence": "TGTTTATTTTG(GT)GGCCTAGC",
    "rightMismatch": "C",
    "reportedMagnesiumMm": 3,
    "loopLength": 39,
    "operator": "Ross",
    "measured": {
      "matched": 73.7,
      "mismatched": 67.7
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 71.8,
        "mismatched": 62.8
      },
      "empiricalWithEnds": {
        "matched": 73.8,
        "mismatched": 65.4
      },
      "santaLucia": {
        "matched": 69.3,
        "mismatched": 60.6
      },
      "rochester": {
        "matched": 64.5,
        "mismatched": 56.2
      }
    }
  },
  {
    "workbookRow": 45,
    "snpName": "SNP19",
    "leftMismatch": "T",
    "coreSequence": "GAGTAGCT(GA)CTAGAATGGGAG",
    "rightMismatch": "A",
    "reportedMagnesiumMm": 3,
    "loopLength": 101,
    "operator": "Hugh",
    "measured": {
      "matched": 70.1,
      "mismatched": 62.9
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 73.6,
        "mismatched": 66
      },
      "empiricalWithEnds": {
        "matched": 76.9,
        "mismatched": 69.5
      },
      "santaLucia": {
        "matched": 72,
        "mismatched": 64.3
      },
      "rochester": {
        "matched": 55.7,
        "mismatched": 47
      }
    }
  },
  {
    "workbookRow": 46,
    "snpName": "SNP19",
    "leftMismatch": "T",
    "coreSequence": "GAGTAGCT(GA)CTAGAATGGGAG",
    "rightMismatch": "A",
    "reportedMagnesiumMm": 2,
    "loopLength": 101,
    "operator": "Ross",
    "measured": {
      "matched": 71.1,
      "mismatched": 64.3
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 72.4,
        "mismatched": 64.8
      },
      "empiricalWithEnds": {
        "matched": 75.8,
        "mismatched": 68.4
      },
      "santaLucia": {
        "matched": 70.7,
        "mismatched": 62.9
      },
      "rochester": {
        "matched": 54.6,
        "mismatched": 46
      }
    }
  },
  {
    "workbookRow": 47,
    "snpName": "SNP20",
    "leftMismatch": "G",
    "coreSequence": "GCTCGGGG(GT)TTAGATAAATAA",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3,
    "loopLength": 53,
    "operator": "Hugh",
    "measured": {
      "matched": 69.1,
      "mismatched": 60.7
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 71.1,
        "mismatched": 63.6
      },
      "empiricalWithEnds": {
        "matched": 74.1,
        "mismatched": 67
      },
      "santaLucia": {
        "matched": 69.7,
        "mismatched": 62.6
      },
      "rochester": {
        "matched": 62.6,
        "mismatched": 55.3
      }
    }
  },
  {
    "workbookRow": 48,
    "snpName": "SNP20",
    "leftMismatch": "G",
    "coreSequence": "GCTCGGGG(GT)TTAGATAAATAA",
    "rightMismatch": "G",
    "reportedMagnesiumMm": 3,
    "loopLength": 53,
    "operator": "Ross",
    "measured": {
      "matched": 72,
      "mismatched": 64.2
    },
    "legacy": {
      "empiricalNoEnds": {
        "matched": 71.1,
        "mismatched": 63.6
      },
      "empiricalWithEnds": {
        "matched": 74.1,
        "mismatched": 67
      },
      "santaLucia": {
        "matched": 69.7,
        "mismatched": 62.6
      },
      "rochester": {
        "matched": 62.6,
        "mismatched": 55.3
      }
    }
  }
]);

export const HUGH_REFIT_ROWS = Object.freeze([
  {
    "workbookRow": 2,
    "experiment": "SNP1  Hugh 4 Mg++",
    "snpName": "SNP1",
    "peak": "Matched",
    "loopLength": 51,
    "recoveredStemTm": 66.26612117624154,
    "originalEmpiricalTm": 79.4
  },
  {
    "workbookRow": 3,
    "experiment": "SNP1  Hugh 4 Mg++",
    "snpName": "SNP1",
    "peak": "Mismatched",
    "loopLength": 51,
    "recoveredStemTm": 61.36767434231083,
    "originalEmpiricalTm": 75.3
  },
  {
    "workbookRow": 4,
    "experiment": "SNP1-7C Hugh 4 Mg++",
    "snpName": "SNP1-7C",
    "peak": "Matched",
    "loopLength": 53,
    "recoveredStemTm": 59.32192152378034,
    "originalEmpiricalTm": 73.5
  },
  {
    "workbookRow": 5,
    "experiment": "SNP1-7C Hugh 4 Mg++",
    "snpName": "SNP1-7C",
    "peak": "Mismatched",
    "loopLength": 53,
    "recoveredStemTm": 51.795039803350235,
    "originalEmpiricalTm": 67.2
  },
  {
    "workbookRow": 6,
    "experiment": "SNP2  Hugh 3 Mg++",
    "snpName": "SNP2",
    "peak": "Matched",
    "loopLength": 24,
    "recoveredStemTm": 50.234299903209006,
    "originalEmpiricalTm": 67.7
  },
  {
    "workbookRow": 7,
    "experiment": "SNP2  Hugh 3 Mg++",
    "snpName": "SNP2",
    "peak": "Mismatched",
    "loopLength": 24,
    "recoveredStemTm": 42.1100466176654,
    "originalEmpiricalTm": 60.9
  },
  {
    "workbookRow": 8,
    "experiment": "SNP3  Hugh 3.5 Mg++",
    "snpName": "SNP3",
    "peak": "Matched",
    "loopLength": 88,
    "recoveredStemTm": 65.72106813475315,
    "originalEmpiricalTm": 77.7
  },
  {
    "workbookRow": 9,
    "experiment": "SNP3  Hugh 3.5 Mg++",
    "snpName": "SNP3",
    "peak": "Mismatched",
    "loopLength": 88,
    "recoveredStemTm": 61.41999286593594,
    "originalEmpiricalTm": 74.1
  },
  {
    "workbookRow": 10,
    "experiment": "SNP4  Hugh 3.5 Mg++",
    "snpName": "SNP4",
    "peak": "Matched",
    "loopLength": 38,
    "recoveredStemTm": 42.28657572549373,
    "originalEmpiricalTm": 60
  },
  {
    "workbookRow": 11,
    "experiment": "SNP4  Hugh 3.5 Mg++",
    "snpName": "SNP4",
    "peak": "Mismatched",
    "loopLength": 38,
    "recoveredStemTm": 35.35706557017713,
    "originalEmpiricalTm": 54.2
  },
  {
    "workbookRow": 12,
    "experiment": "SNP5  Hugh 3.5 Mg++",
    "snpName": "SNP5",
    "peak": "Matched",
    "loopLength": 90,
    "recoveredStemTm": 58.37487834475085,
    "originalEmpiricalTm": 71.5
  },
  {
    "workbookRow": 13,
    "experiment": "SNP5  Hugh 3.5 Mg++",
    "snpName": "SNP5",
    "peak": "Mismatched",
    "loopLength": 90,
    "recoveredStemTm": 51.325893876411534,
    "originalEmpiricalTm": 65.6
  },
  {
    "workbookRow": 14,
    "experiment": "SNP6  Hugh 3.5 Mg++",
    "snpName": "SNP6",
    "peak": "Matched",
    "loopLength": 79,
    "recoveredStemTm": 46.07233241251472,
    "originalEmpiricalTm": 61.5
  },
  {
    "workbookRow": 15,
    "experiment": "SNP6  Hugh 3.5 Mg++",
    "snpName": "SNP6",
    "peak": "Mismatched",
    "loopLength": 79,
    "recoveredStemTm": 40.098616761379716,
    "originalEmpiricalTm": 56.5
  },
  {
    "workbookRow": 16,
    "experiment": "SNP6  Hugh 4 Mg++",
    "snpName": "SNP6",
    "peak": "Matched",
    "loopLength": 79,
    "recoveredStemTm": 46.55022966460552,
    "originalEmpiricalTm": 61.9
  },
  {
    "workbookRow": 17,
    "experiment": "SNP6  Hugh 4 Mg++",
    "snpName": "SNP6",
    "peak": "Mismatched",
    "loopLength": 79,
    "recoveredStemTm": 40.45703970044781,
    "originalEmpiricalTm": 56.8
  },
  {
    "workbookRow": 18,
    "experiment": "SNP7  Hugh 3 Mg++",
    "snpName": "SNP7",
    "peak": "Matched",
    "loopLength": 42,
    "recoveredStemTm": 57.971396385410976,
    "originalEmpiricalTm": 72.9
  },
  {
    "workbookRow": 19,
    "experiment": "SNP7  Hugh 3 Mg++",
    "snpName": "SNP7",
    "peak": "Mismatched",
    "loopLength": 42,
    "recoveredStemTm": 49.01082290870847,
    "originalEmpiricalTm": 65.4
  },
  {
    "workbookRow": 20,
    "experiment": "SNP8  Hugh  3.5 Mg++",
    "snpName": "SNP8",
    "peak": "Matched",
    "loopLength": 74,
    "recoveredStemTm": 59.39482261480003,
    "originalEmpiricalTm": 72.8
  },
  {
    "workbookRow": 21,
    "experiment": "SNP8  Hugh  3.5 Mg++",
    "snpName": "SNP8",
    "peak": "Mismatched",
    "loopLength": 74,
    "recoveredStemTm": 54.13795284180124,
    "originalEmpiricalTm": 68.4
  },
  {
    "workbookRow": 22,
    "experiment": "SNP9  Hugh 3.5 Mg++",
    "snpName": "SNP9",
    "peak": "Matched",
    "loopLength": 71,
    "recoveredStemTm": 64.89737912876362,
    "originalEmpiricalTm": 77.5
  },
  {
    "workbookRow": 23,
    "experiment": "SNP9  Hugh 3.5 Mg++",
    "snpName": "SNP9",
    "peak": "Mismatched",
    "loopLength": 71,
    "recoveredStemTm": 61.07420111203721,
    "originalEmpiricalTm": 74.3
  },
  {
    "workbookRow": 24,
    "experiment": "SNP10  Hugh 3 Mg++",
    "snpName": "SNP10",
    "peak": "Matched",
    "loopLength": 50,
    "recoveredStemTm": 47.45471030198818,
    "originalEmpiricalTm": 63.7
  },
  {
    "workbookRow": 25,
    "experiment": "SNP10  Hugh 3 Mg++",
    "snpName": "SNP10",
    "peak": "Mismatched",
    "loopLength": 50,
    "recoveredStemTm": 34.43201018251386,
    "originalEmpiricalTm": 52.8
  },
  {
    "workbookRow": 26,
    "experiment": "SNP11 Hugh 3 Mg++",
    "snpName": "SNP11",
    "peak": "Matched",
    "loopLength": 73,
    "recoveredStemTm": 65.57042415248793,
    "originalEmpiricalTm": 78
  },
  {
    "workbookRow": 27,
    "experiment": "SNP11 Hugh 3 Mg++",
    "snpName": "SNP11",
    "peak": "Mismatched",
    "loopLength": 73,
    "recoveredStemTm": 57.92406811903512,
    "originalEmpiricalTm": 71.6
  },
  {
    "workbookRow": 28,
    "experiment": "SNP12 Hugh 3 Mg++",
    "snpName": "SNP12",
    "peak": "Matched",
    "loopLength": 83,
    "recoveredStemTm": 58.39326163079377,
    "originalEmpiricalTm": 71.7
  },
  {
    "workbookRow": 29,
    "experiment": "SNP12 Hugh 3 Mg++",
    "snpName": "SNP12",
    "peak": "Mismatched",
    "loopLength": 83,
    "recoveredStemTm": 54.57008361406737,
    "originalEmpiricalTm": 68.5
  },
  {
    "workbookRow": 30,
    "experiment": "SNP13 Hugh 3 Mg++",
    "snpName": "SNP13",
    "peak": "Matched",
    "loopLength": 91,
    "recoveredStemTm": 57.090761421368875,
    "originalEmpiricalTm": 70.4
  },
  {
    "workbookRow": 31,
    "experiment": "SNP13 Hugh 3 Mg++",
    "snpName": "SNP13",
    "peak": "Mismatched",
    "loopLength": 91,
    "recoveredStemTm": 47.29386775350746,
    "originalEmpiricalTm": 62.2
  },
  {
    "workbookRow": 32,
    "experiment": "SNP14 Hugh 3 Mg++",
    "snpName": "SNP14",
    "peak": "Matched",
    "loopLength": 74,
    "recoveredStemTm": 66.56328139616204,
    "originalEmpiricalTm": 78.8
  },
  {
    "workbookRow": 33,
    "experiment": "SNP14 Hugh 3 Mg++",
    "snpName": "SNP14",
    "peak": "Mismatched",
    "loopLength": 74,
    "recoveredStemTm": 62.38168044036753,
    "originalEmpiricalTm": 75.3
  },
  {
    "workbookRow": 34,
    "experiment": "SNP15 Hugh 3 Mg++",
    "snpName": "SNP15",
    "peak": "Matched",
    "loopLength": 37,
    "recoveredStemTm": 53.444514995641256,
    "originalEmpiricalTm": 69.4
  },
  {
    "workbookRow": 35,
    "experiment": "SNP15 Hugh 3 Mg++",
    "snpName": "SNP15",
    "peak": "Mismatched",
    "loopLength": 37,
    "recoveredStemTm": 45.798158962188445,
    "originalEmpiricalTm": 63
  },
  {
    "workbookRow": 36,
    "experiment": "SNP16 Hugh 3 Mg++",
    "snpName": "SNP16",
    "peak": "Matched",
    "loopLength": 82,
    "recoveredStemTm": 60.63025415174971,
    "originalEmpiricalTm": 73.6
  },
  {
    "workbookRow": 37,
    "experiment": "SNP16 Hugh 3 Mg++",
    "snpName": "SNP16",
    "peak": "Mismatched",
    "loopLength": 82,
    "recoveredStemTm": 51.311257735979105,
    "originalEmpiricalTm": 65.8
  },
  {
    "workbookRow": 38,
    "experiment": "SNP17 Hugh 3 Mg++",
    "snpName": "SNP17",
    "peak": "Matched",
    "loopLength": 46,
    "recoveredStemTm": 44.718612444836644,
    "originalEmpiricalTm": 61.6
  },
  {
    "workbookRow": 39,
    "experiment": "SNP17 Hugh 3 Mg++",
    "snpName": "SNP17",
    "peak": "Mismatched",
    "loopLength": 46,
    "recoveredStemTm": 39.46174267183784,
    "originalEmpiricalTm": 57.2
  },
  {
    "workbookRow": 40,
    "experiment": "SNP18 Hugh 3 Mg++",
    "snpName": "SNP18",
    "peak": "Matched",
    "loopLength": 39,
    "recoveredStemTm": 56.45530368803957,
    "originalEmpiricalTm": 71.8
  },
  {
    "workbookRow": 41,
    "experiment": "SNP18 Hugh 3 Mg++",
    "snpName": "SNP18",
    "peak": "Mismatched",
    "loopLength": 39,
    "recoveredStemTm": 45.70261551599656,
    "originalEmpiricalTm": 62.8
  },
  {
    "workbookRow": 42,
    "experiment": "SNP19 Hugh 3 Mg++",
    "snpName": "SNP19",
    "peak": "Matched",
    "loopLength": 101,
    "recoveredStemTm": 61.197953658732224,
    "originalEmpiricalTm": 73.6
  },
  {
    "workbookRow": 43,
    "experiment": "SNP19 Hugh 3 Mg++",
    "snpName": "SNP19",
    "peak": "Mismatched",
    "loopLength": 101,
    "recoveredStemTm": 52.117905869007025,
    "originalEmpiricalTm": 66
  },
  {
    "workbookRow": 44,
    "experiment": "SNP20 Hugh 3 Mg++",
    "snpName": "SNP20",
    "peak": "Matched",
    "loopLength": 53,
    "recoveredStemTm": 56.454538011235535,
    "originalEmpiricalTm": 71.1
  },
  {
    "workbookRow": 45,
    "experiment": "SNP20 Hugh 3 Mg++",
    "snpName": "SNP20",
    "peak": "Mismatched",
    "loopLength": 53,
    "recoveredStemTm": 47.49396453453303,
    "originalEmpiricalTm": 63.6
  }
]);
