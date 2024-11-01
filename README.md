# ARSPA code for testing purposes.

**0. Prerequisites**
- ZoKrates Installation
Refer to the [ZoKrates Getting Started Guide](https://zokrates.github.io/gettingstarted.html) for a one-line installation command.

- Setting up ZoKrates Path
To ensure the ZoKrates binary is accessible from any directory, add it to your system's PATH. Run the following command to add the path to your shell’s configuration file (e.g., .bashrc, .zshrc), making the change permanent:
  ```
  echo 'export PATH=$PATH:$HOME/.zokrates/bin' >> ~/.zshrc  # or ~/.bashrc
  source ~/.zshrc  # apply the change
  ```
 ### Remix IDE Setup

To use Remix IDE for smart contract development, refer to the official [Remix IDE Documentation](https://remix-ide.readthedocs.io/en/latest/) for setup instructions and configuration tips.

You can access Remix IDE directly in your browser by visiting [Remix IDE](https://remix.ethereum.org/). This IDE is highly recommended for writing, testing, and deploying smart contracts with ease.

#### Basic Steps:
1. Open [Remix IDE](https://remix.ethereum.org/) in your browser.
2. Select the appropriate environment for your development needs (e.g., JavaScript VM, Injected Web3).
3. Start coding, compiling, and testing your smart contracts directly within the IDE.

For additional plugins and advanced configurations, refer to the [Remix Plugin Manager](https://remix-ide.readthedocs.io/en/latest/plugin_manager.html).

---

이렇게 작성하면, Remix IDE를 사용하는 방법과 관련 링크를 포함하여 사용자에게 필요한 설정 정보를 명확하게 전달할 수 있습니다.
 
**1. Smart contract**
- The smart contract can be deployed to the blockchain. Before deployment, the zk-SNARK verification keys inside the smart contract should be changed to the zk-SNARK keys generated by the HLDS programs.

**2. zk-SNARK programs in Zokrates**
- Generate zk-SNARK keypair: (the verification keys should be embedded to the smart contract)
```
zokrates compile -i <name>.zok
zokrates set -s gm17
```
- Generate proof:
```
zokrates compute-witness -a <inputs of name.zok>
zokrates generate-proof -s gm17
zokrates print-proof --format remix
```
- The output is the proof string $\pi$ and public input/output ($\vec{v},\vec{o}$). They can be passed directly to the smart contract functions $\mathsf{verifyTx}$.

**3. zkSNARK programs' inputs generation**
- The root and leaves of the Merkle tree can be computed using Zokrates programs in the input folders (input7, input10, input13, input17)
  ```
  zokrates compile -i <name>.zok
  zokrates compute-witness --verbose >> output.txt
  ```
- For example, to compute a leaf ($\mathsf{leaf = H(H(r),H(cID))}$), we can run functions $c<number>$:
  ```
  zokrates compile -i c1.zok
  zokrates compute-witness --verbose >> output.txt
  ```
- The leaf is written to output.txt:
  ```
  Computing witness...

  Witness: 
  ["12068829434966858217616654989280481420045420412202930300116223112089659876982"]

  Witness file written to 'witness'
  ```
- Similarly, we can compute the root by using functions $r<number>$
  ```
  zokrates compile -i r1.zok
  zokrates compute-witness --verbose >> output.txt
  ```
- The root is written to output.txt:
  ```
  Computing witness...

  Witness: 
  ["17088704704177241315644229239628817314264008424491911849361778587561865360994"]

  Witness file written to 'witness'
  ```
**4. generate zkSNARK proof**
- The main zokrates programs are in folders: depth7 ($2^7$ users), depth10 ($2^{10}$ users), depth13 ($2^{13}$ users), and depth17 ($2^{17}$ users).
- In each folder, there are 13 programs corresponding to the number of candidates (from 1 candidate to 13 candidates each reviewer).
- For example, to generate a proof for a system of $2^{17}$ users and each user reviews 13 candidates, we can use /arspa/depth17/c13/c13_d17.zok:
  ```
  zokrates compile -i c13_d17.zok
  zokrates setup -s gm17
  zokrates compute witness -a <inputs of c13_d17.zok generated from the above steps>
  zokrates generate-proof -s gm17
  zokrates print-proof --format remix >> proof_c13d17.txt
  ```
- The verification key in the smart contract should be updated using the newly generated verification key.
- The proof string in proof_c13d17.txt can be used to call the function verifyTx in the smart contract.

**5. Smart contract**
- The smart contract can be deployed in the main net or testnet.
- We can test the verifyTx() function by using the proof string from the above steps.
- NOTE: candidate IDs must be in the range of 0 to 12.
