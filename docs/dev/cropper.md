## Algorithm 3: AFFINITY CROPPER
Come from BOLTZ-2 Algorithm 3: Affinity Cropper
**Input:**
- Token list `tokens`
- Minimum distance to the ligand `min_dist_protein_to_ligand`
- `max_tokens = 256`, `max_protein_tokens = 200`, `neighborhood_size = 10`

```text
// Start with all ligand tokens
cropped_tokens ← tokens[mol_type = ligand]

// Min pooling over the residues' atoms
min_dist_res_to_ligand ← MinPooling(min_dist_protein_to_ligand[protein])
res_idx_sorted ← argsort(min_dist_res_to_ligand)

// Add protein tokens around pocket residues
for res_idx in res_idx_sorted do:
    Let res_tokens be the entries with res_idx
    Let chain_id be the asym_id of the current residue
    Let chain_tokens be protein tokens with asym_id = chain_id

    // Initialize residue window
    min_idx = max_idx = res_idx
    while len(res_tokens) < neighborhood_size do:
        min_idx = min_idx - 1
        max_idx = max_idx + 1
        res_tokens ← all tokens in chain_tokens with res_idx ∈ [min_idx, max_idx]

    Let new_tokens be the entries in res_tokens not in cropped_tokens

    // Check token limits
    if (cropped_tokens ∪ new_tokens > max_tokens) or
       ((cropped_tokens ∪ new_tokens) ∩ protein_tokens > max_protein_tokens) then:
        break

    cropped_tokens ← cropped_tokens ∪ new_tokens
end for

Output: cropped_tokens
```

---

### 格式说明：

- 使用 **代码块**（```text）包裹伪代码，保持等宽字体和缩进。
- 注释使用 `//` 开头，与原图一致。
- 使用 `←` 表示赋值（可替换为 `=` 若平台不支持 Unicode）。
- 集合操作符如 `∪`（并集）、`∈`（属于）、`∩`（交集）保留 Unicode 符号，兼容主流 Markdown 渲染器。
- 缩进使用 4 个空格，清晰表达控制结构层级。
- 条件判断中的多行 `if` 使用换行+缩进，提高可读性。

---

✅ 此版本可直接复制粘贴到支持 Markdown 的编辑器中使用。如需更简洁的纯文本版或转换为 HTML/Word/PDF，也可告知我为您生成。
