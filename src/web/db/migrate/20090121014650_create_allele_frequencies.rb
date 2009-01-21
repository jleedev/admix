class CreateAlleleFrequencies < ActiveRecord::Migration
  def self.up
    create_table :allele_frequencies do |t|
      t.float :freq
      t.integer :population_number
      t.belongs_to :marker

      t.timestamps
    end
  end

  def self.down
    drop_table :allele_frequencies
  end
end
