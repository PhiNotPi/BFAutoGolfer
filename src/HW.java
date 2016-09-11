import java.util.*;

public class HW {

  // public static enum Command {
  // L, R, INC, DEC, OB, CB, IN, OUT
  // }

  // static final Command[] commandlist = { Command.L, Command.R, Command.INC,
  // Command.DEC, Command.OB, Command.CB, Command.IN, Command.OUT };

  static final int L = 1;
  static final int R = 2;
  static final int INC = 3;
  static final int DEC = 4;
  static final int OB = 5;
  static final int CB = 6;
  static final int IN = 7;
  static final int OUT = 8;

  public static class CommandPair {
    int type;
    int run;

    public CommandPair(int com, int run) {
      this.type = com;
      this.run = run;
    }

    public String toString() {
      String res = "";
      String typeStr = commandToString(type);
      for (int i = 0; i < run; i++) {
        res += typeStr;
      }
      return res;
    }
  }

  public static String commandToString(int type) {
    switch (type) {
    case L:
      return "<";
    case R:
      return ">";
    case INC:
      return "+";
    case DEC:
      return "-";
    case OB:
      return "[";
    case CB:
      return "]";
    case IN:
      return ",";
    case OUT:
      return ".";
    }
    return "?";
  }

  public static class CommandPairList {
    private ArrayList<CommandPair> coms = new ArrayList<CommandPair>();
    int lengthSum = 0;

    public void add(CommandPair com) {
      coms.add(com);
      lengthSum += com.run;
    }

    public void add(int loc, CommandPair com) {
      coms.add(loc, com);
      lengthSum += com.run;
    }

    public void addAll(CommandPairList data) {
      this.addAll(this.coms.size(), data);
    }

    public void addAll(int loc, CommandPairList data) {
      coms.addAll(loc, data.coms);
      lengthSum += data.lengthSum;
    }

    public void addAll(int[] data) {
      this.addAll(this.coms.size(), data);
    }

    public void addAll(int loc, int[] data) {
      CommandPairList toAdd = new CommandPairList();
      for (int i = 0; i < data.length; i++) {
        if (toAdd.coms.size() == 0
            || toAdd.get(toAdd.coms.size() - 1).type != data[i]) {
          toAdd.add(new CommandPair(data[i], 0));
        }
        toAdd.get(toAdd.coms.size() - 1).run++;
        toAdd.lengthSum++;
      }
      this.addAll(loc, toAdd);
    }

    public CommandPair get(int loc) {
      return coms.get(loc);
    }

    public String toString() {
      String res = lengthSum + ": ";
      for (CommandPair com : coms) {
        res += com;
      }
      return res;
    }

    public int[] toCommandArray() {
      int[] res = new int[lengthSum];
      int i = 0;
      for (CommandPair com : coms) {
        for (int j = 0; j < com.run; j++, i++) {
          res[i] = com.type;
        }
      }
      return res;
    }
  }

  public static final byte[] targetString = stringToByteArray("PhiNotPi");

  public static void main(String[] args) {
    // TODO
    // byte[] target = stringToByteArray("Hello, World!");
    // CommandPairList prog = simplePathfinder(new byte[15], 5, target, 0,
    // Integer.MAX_VALUE);
    // System.out.println(eval(prog.toCommandArray()));
    // cycleLengths(15, 30);
    // System.out.println(isValidRecRel(new int[] { DEC, DEC, L, DEC, L, L, INC,
    // OB, INC, OB, L, INC, R, DEC, DEC, DEC, R, DEC, R, DEC, R, DEC, L, L, L,
    // CB, R, CB }));
    // printNestedLoops("abcdefghijklmnop");
    iterRecRel(28);
  }

  public static byte[] stringToByteArray(String s) {
    byte[] res = new byte[s.length()];
    for (int i = 0; i < s.length(); i++) {
      res[i] = (byte) s.charAt(i);
    }
    return res;
  }

  static CommandPairList simplePathfinder(final byte[] tape, final int loc,
      final byte[] reqout, final int curchar, int best) {
    CommandPairList bestSuffix = null;
    for (Integer i = loc; i != null; i = staggerNext(i, loc, tape.length)) {
      CommandPairList curSuffix = new CommandPairList();
      if (i < loc) {
        curSuffix.add(new CommandPair(L, loc - i));
      } else if (i > loc) {
        curSuffix.add(new CommandPair(R, i - loc));
      }
      int dif = (byte) (reqout[curchar] - tape[i]);
      int absdif = dif;
      if (absdif < 0) {
        absdif = -absdif;
      }
      int absdest = reqout[curchar];
      if (absdest < 0) {
        absdest = -absdest;
      }
      if (absdif <= absdest + 1) {
        if (dif > 0) {
          curSuffix.add(new CommandPair(INC, absdif));
        } else if (dif < 0) {
          curSuffix.add(new CommandPair(DEC, absdif));
        }
      } else {
        curSuffix.add(new CommandPair(IN, 1));
        if (reqout[curchar] > 0) {
          curSuffix.add(new CommandPair(INC, absdest));
        } else if (reqout[curchar] < 0) {
          curSuffix.add(new CommandPair(DEC, absdest));
        }
      }
      curSuffix.add(new CommandPair(OUT, 1));
      if (curSuffix.lengthSum >= best) {
        continue;
      }
      if (curchar < reqout.length - 1) {
        byte[] childtape = tape.clone();
        childtape[i] = reqout[curchar];
        CommandPairList childSuf = simplePathfinder(childtape, i, reqout,
            curchar + 1, best - curSuffix.lengthSum);
        if (childSuf == null) {
          continue;
        }
        curSuffix.addAll(childSuf);
      }
      // if (curchar == 0) {
      // System.out.println(curSuffix);
      // }
      best = curSuffix.lengthSum;
      bestSuffix = curSuffix;
    }
    return bestSuffix;
  }

  public static Integer staggerNext(int cur, int pivot, int limit) {
    Integer res = null;
    if (cur > pivot) {
      res = (pivot << 1) - cur;
      if (res < 0) {
        res = cur + 1;
        if (res >= limit) {
          return null;
        }
      }
    } else {
      res = (pivot << 1) - cur + 1;
      if (res >= limit) {
        res = cur - 1;
        if (res < 0) {
          return null;
        }
      }
    }
    return res;
  }

  public static class TapeState {
    ArrayList<Byte> tape;
    ArrayList<Character> output;
    int dp;
    int time;

    public TapeState() {
      tape = new ArrayList<Byte>();
      output = new ArrayList<Character>();
      dp = 0;
      time = 0;
    }

    public String toString() {
      return output + "" + tape + dp;
    }

    public byte[] tapeArray() {
      byte[] res = new byte[tape.size()];
      for (int i = 0; i < tape.size(); i++) {
        res[i] = tape.get(i);
      }
      return res;
    }

    public byte[] toByteArray() {
      // byte[] res = new byte[8+output.size()+tape.size()];
      byte[] res = new byte[4 + tape.size()];
      int j = 0;
      res[j++] = (byte) ((dp >>> 0) & 0xff);
      res[j++] = (byte) ((dp >>> 8) & 0xff);
      res[j++] = (byte) ((dp >>> 16) & 0xff);
      res[j++] = (byte) ((dp >>> 24) & 0xff);
      // res[j++] = (byte) ((output.size() >>> 0) & 0xff);
      // res[j++] = (byte) ((output.size() >>> 8) & 0xff);
      // res[j++] = (byte) ((output.size() >>> 16) & 0xff);
      // res[j++] = (byte) ((output.size() >>> 24) & 0xff);
      // for(int i = 0; i < output.size(); i++){
      // res[j++] = (byte)(output.get(i).charValue());
      // }
      for (int i = 0; i < tape.size(); i++) {
        res[j++] = tape.get(i);
      }
      return res;
    }
  }

  public static TapeState eval(int[] coms) {
    TapeState tapestate = new TapeState();
    ArrayList<Byte> tape = tapestate.tape;
    ArrayList<Character> output = tapestate.output;
    Map<Integer, Integer> match = new HashMap<Integer, Integer>();
    Stack<Integer> bracestack = new Stack<Integer>();
    for (int i = 0; i < coms.length; i++) {
      if (coms[i] == OB) {
        bracestack.push(i);
      } else if (coms[i] == CB) {
        if (bracestack.isEmpty()) {
          return null;
        }
        int matchindex = bracestack.pop();
        match.put(i, matchindex);
        match.put(matchindex, i);
      }
    }
    if (!bracestack.isEmpty()) {
      return null;
    }
    int dp = 0;
    tape.add(0, (byte) 0);
    for (int ip = 0; ip < coms.length; ip++, tapestate.time++) {
      if (tapestate.time > 100000 || tape.size() > 300) {
        return null;
      }
      switch (coms[ip]) {
      case L:
        if (dp == 0) {
          tape.add(0, (byte) 0);
        } else {
          dp--;
        }
        break;
      case R:
        dp++;
        if (dp == tape.size()) {
          tape.add((byte) 0);
        }
        break;
      case INC:
        tape.set(dp, (byte) (tape.get(dp) + 1));
        break;
      case DEC:
        tape.set(dp, (byte) (tape.get(dp) - 1));
        break;
      case OB:
        if (tape.get(dp) == 0) {
          ip = match.get(ip);
        }
        break;
      case CB:
        if (tape.get(dp) != 0) {
          ip = match.get(ip);
        }
        break;
      case IN:
        tape.set(dp, (byte) 0);
        break;
      case OUT:
        output.add((char) (tape.get(dp).byteValue()));
        break;
      }
    }
    tapestate.dp = dp;
    return tapestate;
  }

  public static CommandPairList cycleLengths(int min, int max) {
    int totalToBeat = Integer.MAX_VALUE;
    CommandPairList bestResult = null;
    BloomFilter<TapeState> bf = new BloomFilter<TapeState>(0.01, 100000);
    for (int i = min; i <= max; i++) {
      CommandPairList lenRes = cycleWithLength(i, totalToBeat, bf);
      if (lenRes != null) {
        totalToBeat = lenRes.lengthSum;
        bestResult = lenRes;
      }
    }
    return bestResult;
  }

  public static CommandPairList cycleWithLength(int length, int totalToBeat,
      BloomFilter<TapeState> bf) {
    int[] commands = new int[length];
    for (int i = 0; i < length - 1; i++) {
      commands[i] = 1;
    }
    commands[length - 1] = 6;
    int bfcounter = 0;
    int dupecounter = 0;
    long invalidcounter = 0;
    long totalcounter = 0;
    long expectedtotal = (long) Math.pow(6, length - 1);
    System.out.println(expectedtotal);
    int percentDone = 0;
    int infloopcounter = 0;
    int bestSuffixLength = totalToBeat - length;
    CommandPairList bestResult = null;
    while (true) {
      if (isValidRecRel(commands)) {
        TapeState evalRes = eval(commands);
        if (evalRes == null) {
          infloopcounter++;
          // System.out.print("{");
          // for (int i = 0; i < length; i++) {
          // System.out.print(commandToString(prog[i]));
          // }
          // System.out.println("}");
        } else {
          if (bf.contains(evalRes.toByteArray())) {
            dupecounter++;
          } else {
            CommandPairList suffix = simplePathfinder(evalRes.tapeArray(),
                evalRes.dp, targetString, 0, bestSuffixLength);
            if (suffix != null) {
              System.out.print("{");
              for (int i = 0; i < length; i++) {
                System.out.print(commandToString(commands[i]));
              }
              System.out.print("}");
              System.out.println(evalRes);
              System.out.println(suffix);
              bestSuffixLength = suffix.lengthSum;
              bestResult = suffix;
              bestResult.addAll(0, commands);
            }
            bfcounter++;
            bf.add(evalRes.toByteArray());
          }
        }
      } else {
        invalidcounter++;
      }
      int i = 0;
      while (commands[i] == 6) {
        commands[i] = 1;
        i++;
        if (i == length) {
          System.out.println("len:" + length + " unique:" + bfcounter
              + " dupe:" + dupecounter + " infloop:" + infloopcounter
              + " invalid:" + invalidcounter);
          System.out.println("BEST: " + bestResult);
          return bestResult;
        }
      }
      commands[i]++;
      totalcounter++;
      if (totalcounter * 100.0 / expectedtotal > percentDone + 5) {
        percentDone += 5;
        System.out.println(percentDone + "%");
      }
    }
  }

  // public static boolean areValidCommands(int[] commands) {
  // int bracecount = 0;
  // int bracepaircount = 0;
  // if (commands.length == 0) {
  // return false;
  // }
  // if (commands[0] != INC && commands[0] != DEC) {
  // return false;
  // }
  // // requirements for usefulness for HW
  // if (commands[commands.length - 1] != CB) {
  // return false;
  // }
  //
  // int sign = 1;
  // int delta = 0;
  // int i = 0;
  // while (i < commands.length) {
  // if (commands[i] == INC) {
  // do {
  // i++;
  // } while (commands[i] == INC);
  // if (commands[i] == DEC) {
  // return false;
  // } else if (commands[i] == R) {
  // do {
  // i++;
  // delta++;
  // } while (commands[i] == R);
  // if (commands[i] == L) {
  // return false;
  // }
  // if (sign == 0) {
  // sign = 1;
  // }
  // if (sign == -1) {
  // if (commands[i] != OB && commands[i] != CB) {
  // return false;
  // }
  // }
  // } else if (commands[i] == L) {
  // do {
  // i++;
  // delta--;
  // } while (commands[i] == L);
  // if (commands[i] == R) {
  // return false;
  // }
  // if (sign == 0) {
  // sign = -1;
  // }
  // if (sign == 1) {
  // if (commands[i] != OB && commands[i] != CB) {
  // return false;
  // }
  // }
  // }
  // } else if (commands[i] == DEC) {
  // do {
  // i++;
  // } while (commands[i] == DEC);
  // if (commands[i] == INC) {
  // return false;
  // } else if (commands[i] == R) {
  // do {
  // i++;
  // delta++;
  // } while (commands[i] == R);
  // if (commands[i] == L) {
  // return false;
  // }
  // if (sign == 0) {
  // sign = 1;
  // }
  // if (sign == -1) {
  // if (commands[i] != OB && commands[i] != CB) {
  // return false;
  // }
  // }
  // } else if (commands[i] == L) {
  // do {
  // i++;
  // delta--;
  // } while (commands[i] == L);
  // if (commands[i] == R) {
  // return false;
  // }
  // if (sign == 0) {
  // sign = -1;
  // }
  // if (sign == 1) {
  // if (commands[i] != OB && commands[i] != CB) {
  // return false;
  // }
  // }
  // } else {
  // return false;
  // }
  // } else if (commands[i] == OB) {
  // if (sign * delta < 0) {
  // return false;
  // }
  // sign = 0;
  // delta = 0;
  // do {
  // bracecount++;
  // bracepaircount++;
  // i++;
  // } while (commands[i] == OB);
  // if (commands[i] == CB || bracepaircount > 2) {
  // return false;
  // }
  // } else if (commands[i] == CB) {
  // if (sign * delta < 0) {
  // return false;
  // }
  // sign = 0;
  // delta = 0;
  // bracecount--;
  // i++;
  // if (bracecount < 0) {
  // return false;
  // }
  // } else if (commands[i] == R) {
  // do {
  // i++;
  // delta++;
  // } while (commands[i] == R);
  // if (commands[i] == L) {
  // return false;
  // }
  // } else if (commands[i] == L) {
  // do {
  // i++;
  // delta--;
  // } while (commands[i] == L);
  // if (commands[i] == R) {
  // return false;
  // }
  // }
  // }
  // if (bracecount != 0) {
  // return false;
  // }
  // return true;
  // }

  public static boolean isValidRecRel(int[] commands) {
    if (commands.length == 0) {
      return false;
    }
    if (commands[0] != INC && commands[0] != DEC) {
      return false;
    }
    // requirements for usefulness for HW
    if (commands[commands.length - 1] != CB) {
      return false;
    }

    int i = 0;
    while (true) {
      if (commands[i] == INC || commands[i] == DEC) {
        int type = commands[i];
        do {
          i++;
        } while (commands[i] == type);
        if (commands[i] == L) {
          do {
            i++;
          } while (commands[i] == L);
          continue;
        } else if (commands[i] == OB) {
          i++;
          break;
        } else {
          return false;
        }
      } else {
        return false;
      }
    }
    // System.out.println("a"+i);
    if (commands[i] == INC || commands[i] == DEC) {
      int type = commands[i];
      do {
        i++;
      } while (commands[i] == type);
    }
    if (commands[i] != OB) {
      return false;
    }
    i++;

    // <+>-
    if (commands[i] != L) {
      return false;
    }
    i++;
    if (commands[i] == INC || commands[i] == DEC) {
      int type = commands[i];
      do {
        i++;
      } while (commands[i] == type);
    } else {
      return false;
    }
    if (commands[i] != R) {
      return false;
    }
    i++;
    if (commands[i] == INC || commands[i] == DEC) {
      int type = commands[i];
      do {
        i++;
      } while (commands[i] == type);
    } else {
      return false;
    }

    // System.out.println("b" + i);
    int pastrun = 0;
    int delta = 0;
    while (true) {
      if (commands[i] == R) {
        pastrun = 0;
        do {
          i++;
          delta++;
        } while (commands[i] == R);
        if (commands[i] == INC || commands[i] == DEC) {
          int type = commands[i];
          do {
            i++;
            pastrun++;
          } while (commands[i] == type);
        } else {
          return false;
        }
      } else if (commands[i] == L) {
        if (pastrun == 0) {
          return false;
        }
        do {
          i++;
          delta--;
        } while (commands[i] == L);
        if (delta != 0) {
          return false;
        }
        break;
      } else {
        return false;
      }
    }

    if (commands[i] != CB || i >= commands.length - 1) {
      return false;
    }
    i++;

    // {}>{k}]
    if (commands[i] == INC || commands[i] == DEC) {
      int type = commands[i];
      do {
        i++;
      } while (commands[i] == type);
    }
    if (commands[i] != R) {
      return false;
    }
    i++;
    if (commands[i] == INC || commands[i] == DEC) {
      int type = commands[i];
      do {
        i++;
      } while (commands[i] == type);
    }
    if (commands[i] != CB) {
      return false;
    }
    i++;
    if (i != commands.length) {
      return false;
    }

    return true;
  }

  // public static boolean isValidPair(int i, int j) {
  // if ((i == INC && j == DEC) || (i == DEC && j == INC)) {
  // return false;
  // }
  // if ((i == L && j == R) || (i == R && j == L)) {
  // return false;
  // }
  // if (i == OB && j == CB) {
  // return false;
  // }
  // return true;
  // }

  public static void printNestedLoops(String s) {
    for (int i = 0; i < s.length(); i++) {
      char curchar = s.charAt(i);
      if (i == 0) {
        System.out.println("for(int " + curchar + " = 0; " + curchar
            + " < length; " + curchar + "++){");
      } else {
        char prevchar = s.charAt(i - 1);
        System.out.println("for(int " + curchar + " = 0; " + curchar + " < "
            + prevchar + "; " + curchar + "++){");
      }
    }
    for (int i = 0; i < s.length(); i++) {
      System.out.print("}");
    }
    System.out.println();
  }

  public static void iterRecRel(int maxlength) {

    BloomFilter<TapeState> bf = new BloomFilter<TapeState>(0.01, 100000);
    long dupecounter = 0;
    long bfcounter = 0;
    long infloopcounter = 0;
    CommandPairList bestResult = null;
    int bestSuffixLength = Integer.MAX_VALUE;

    for (int a = 0; a < maxlength; a++) {
      for (int b = 0; b < a; b++) {
        for (int c = 0; c < b; c++) {
          for (int d = c - 1; d < c; d++) {
            for (int e = d - 1; e < d; e++) {
              for (int f = e - 1; f < e; f++) {
                for (int g = 0; g < f; g++) {
                  for (int h = 0; h < g; h++) {
                    for (int i = 0; i < h; i++) {
                      for (int j = i % 2; j < i - 1; j += 2) {
                        for (int k = 0; k < j - 1; k++) {
                          for (int l = k - 1; l < k; l++) {
                            for (int m = 0; m < l; m++) {
                              for (int n = 0; n < m - 1; n++) {
                                for (int o = 0; o < n; o++) {
                                  for (int p = 1; p < o; p++) {
                                    int[] commands = new int[a + 1];
                                    // --<-<<+[+[<+>--->->->-<<<]>]
                                    // <<<[[<>>>><<<]>]
                                    commands[a] = CB;
                                    commands[b] = R;
                                    commands[c] = CB;
                                    commands[d] = L;
                                    commands[e] = L;
                                    commands[f] = L;
                                    commands[g] = R;
                                    commands[h] = R;
                                    commands[i] = R;
                                    commands[j] = R;
                                    commands[k] = L;
                                    commands[l] = OB;
                                    commands[m] = OB;
                                    commands[n] = L;
                                    commands[o] = L;
                                    commands[p] = L;
                                    int runcount = 0;
                                    boolean inrun = false;
                                    for (int ip = 0; ip < commands.length; ip++) {
                                      if (commands[ip] == 0) {
                                        if (!inrun) {
                                          inrun = true;
                                          runcount++;
                                        }
                                      } else {
                                        inrun = false;
                                      }
                                      // System.out
                                      // .print(commandToString(commands[ip]));
                                    }
                                    // System.out.println();

                                    int numoptions = 1 << runcount;
                                    for (int curoption = 0; curoption < numoptions; curoption++) {
                                      int[] prog = commands.clone();

                                      runcount = 0;
                                      inrun = false;
                                      for (int ip = 0; ip < prog.length; ip++) {
                                        if (prog[ip] == 0) {
                                          if (!inrun) {
                                            inrun = true;
                                            runcount++;
                                          }
                                          if (((curoption >>> (runcount - 1)) & 1) == 1) {
                                            prog[ip] = INC;
                                          } else {
                                            prog[ip] = DEC;
                                          }
                                        } else {
                                          inrun = false;
                                        }
                                        // System.out
                                        // .print(commandToString(prog[ip]));
                                      }
                                      // System.out.println();

                                      TapeState evalRes = eval(prog);
                                      if (evalRes == null) {
                                        infloopcounter++;
                                      } else {
                                        if (bf.contains(evalRes.toByteArray())) {
                                          dupecounter++;
                                        } else {
                                          CommandPairList suffix = simplePathfinder(
                                              evalRes.tapeArray(), evalRes.dp,
                                              targetString, 0, bestSuffixLength
                                                  - prog.length);
                                          if (suffix != null) {
                                            System.out.print("{");
                                            for (int ip = 0; ip < prog.length; ip++) {
                                              System.out
                                                  .print(commandToString(prog[ip]));
                                            }
                                            System.out.print("}");
                                            System.out.println(evalRes);
                                            bestResult = suffix;
                                            bestResult.addAll(0, prog);
                                            System.out.println(bestResult);
                                            System.out.println("len:"
                                                + prog.length + " unique:"
                                                + bfcounter + " dupe:"
                                                + dupecounter + " inf:"
                                                + infloopcounter);
                                            bestSuffixLength = bestResult.lengthSum;
                                          }
                                          bfcounter++;
                                          bf.add(evalRes.toByteArray());
                                        }
                                      }

                                    }

                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

}
